#include "filter.h"
#include "../solve/solver.h"

#define BUCKET_INDEX(x,y,z,w,h,d) ((x) * (h) * (d) + (y) * (d) + (z))

#ifdef SUPPORT_GPU
__global__ void filter_gaussian_kernel(const Precision* xyz,
                                  const Precision* val,
                                  const ID* ar1,
                                  const ID* ar2,
                                        Precision* res,
                                  const Precision rad,
                                  const Precision sig,
                                  const ID n_nodes,
                                  const ID blocks_x,
                                  const ID blocks_y,
                                  const ID blocks_z) {
    int n_id = blockIdx.x * blockDim.x + threadIdx.x;

    if(n_id >= n_nodes) return;

    // extract location
    Precision x = xyz[3 * n_id + 0];
    Precision y = xyz[3 * n_id + 1];
    Precision z = xyz[3 * n_id + 2];

    // compute the block location
    ID bucket_x = x / rad;
    ID bucket_y = y / rad;
    ID bucket_z = z / rad;

    Precision weight = 0;
    Precision value = 0;
    for(int bx = max(bucket_x - 1, 0); bx < min(bucket_x + 2, blocks_x); bx++){
        for(int by = max(bucket_y - 1, 0); by < min(bucket_y + 2, blocks_y); by++){
            for(int bz = max(bucket_z - 1, 0); bz < min(bucket_z + 2, blocks_z); bz++){
                ID id_start = ar1[BUCKET_INDEX(bx,by,bz, blocks_x, blocks_y, blocks_z)];
                ID id_end   = ar1[BUCKET_INDEX(bx,by,bz, blocks_x, blocks_y, blocks_z) + 1];
                for(ID i = id_start; i < id_end; i++){
                    ID n_id_2 = ar2[i];

                    Precision x_2 = xyz[3 * n_id_2 + 0];
                    Precision y_2 = xyz[3 * n_id_2 + 1];
                    Precision z_2 = xyz[3 * n_id_2 + 2];

                    Precision dist_squared = (x-x_2) * (x-x_2) + (y-y_2) * (y-y_2) + (z-z_2) * (z-z_2);

                    if(dist_squared < rad * rad){
                        Precision wgt = exp(- dist_squared / (2 * sig * sig));
                        weight += wgt;
                        value += wgt * val[n_id_2];
                    }
                }
            }
        }
    }
    res[n_id] = value / weight;
}
#endif

DynamicVector fem::math::filter_gaussian(DynamicMatrix &coords, DynamicVector &values, Precision radius, Precision sigma, solver::SolverDevice device){
    // dont do anything if radius or sigma is 0
    if(radius == 0 || sigma == 0){
        return values;
    }

    // step 1, copy coordinates to perform operations
    DynamicMatrix xyz = DynamicMatrix(coords);

    // step 2, compute bounds
    Precision min_x = xyz.col(0).minCoeff();
    Precision min_y = xyz.col(1).minCoeff();
    Precision min_z = xyz.col(2).minCoeff();
    Precision max_x = xyz.col(0).maxCoeff();
    Precision max_y = xyz.col(1).maxCoeff();
    Precision max_z = xyz.col(2).maxCoeff();

    // step 3, adjust xyz so that min_x and so on are at 0,0,0
    xyz.col(0) = xyz.col(0).array() - min_x;
    xyz.col(1) = xyz.col(1).array() - min_y;
    xyz.col(2) = xyz.col(2).array() - min_z;
    max_x -= min_x;
    max_y -= min_y;
    max_z -= min_z;

    // adjust radius based on a minimal size constraint
    Precision volume = max_x * max_y * max_z;
    Precision vol_pp = volume / xyz.rows();
    Precision scale  = std::cbrt(vol_pp);
    radius = std::max(scale / 3, radius);

    // step 4, compute grid size
    size_t grid_x  = std::ceil(max_x / radius);
    size_t grid_y  = std::ceil(max_y / radius);
    size_t grid_z  = std::ceil(max_z / radius);
    size_t buckets = grid_x * grid_y * grid_z;

    logging::info(true, "Filter radius   : ", radius);
    logging::info(true, "Filter grid size: ", grid_x,", ",grid_y, ", ", grid_z);
    logging::info(true, "Filter buckets  : ", buckets);

    // step 5, allocate arrays for the grid
    // for this purpose we use two arrays. the first one maps bucket index to offset at which the first node is
    // the second array contains the node ids for the given bucket. for example if there are 2 buckets, the first one
    // with ids 1,3,5 and the second one with ids 2 and 4:
    // array_1 = [0, 3, 5]
    // array_2 = [1, 3, 5, 2, 4]
    IndexVector array_1 = IndexVector::Zero(buckets + 1);
    IndexVector array_2 = IndexVector::Zero(coords.rows());

    // step 6, compute how many times a bucket is used
    IndexVector bucket_node_counts = IndexVector::Zero(buckets);
    for(int m = 0; m < xyz.rows(); m++){
        ID bucket_x = xyz(m, 0) / radius;
        ID bucket_y = xyz(m, 1) / radius;
        ID bucket_z = xyz(m, 2) / radius;
        ID bucket_idx = BUCKET_INDEX(bucket_x, bucket_y, bucket_z, grid_x, grid_y, grid_z);
        bucket_node_counts(bucket_idx) ++;
    }

    // step 7, create a temporary large matrix which holds the indices
    IndexMatrix bucket_temp = IndexMatrix::Constant(buckets, bucket_node_counts.maxCoeff(), -1);
    for(int m = 0; m < xyz.rows(); m++){
        ID bucket_x = xyz(m, 0) / radius;
        ID bucket_y = xyz(m, 1) / radius;
        ID bucket_z = xyz(m, 2) / radius;
        ID bucket_idx = BUCKET_INDEX(bucket_x, bucket_y, bucket_z, grid_x, grid_y, grid_z);
        bucket_node_counts(bucket_idx) --;
        bucket_temp(bucket_idx, bucket_node_counts(bucket_idx)) = m;
    }

    // step 8, convert to the previously declared arrays
    for(int m = 0; m < buckets; m++){
        int c = 0;
        for(c = 0; c < bucket_temp.cols(); c++){
            if(bucket_temp(m,c) < 0) break;
            else array_2(array_1(m) + c) = bucket_temp(m,c);
        }
        array_1(m+1) = array_1(m) + c;
    }
    // perform operations, either on gpu if gpu is supported, else the cpu
#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#endif

#ifdef SUPPORT_GPU
    if(device == solver::GPU){
        // create CudaArray objects
        cuda::CudaArray<Precision> d_xyz(coords.rows() * 3);
        cuda::CudaArray<Precision> d_val(values.size());
        cuda::CudaArray<ID>        d_ar1(array_1.size());
        cuda::CudaArray<ID>        d_ar2(array_2.size());
        cuda::CudaArray<Precision> d_res(values.size());

        // upload data to GPU
        d_xyz.upload(xyz.data());
        d_val.upload(values.data());
        d_ar1.upload(array_1.data());
        d_ar2.upload(array_2.data());

        // define block size and grid size
        dim3 blockSize(256);
        dim3 gridSize((values.size() + blockSize.x - 1) / blockSize.x);

        // call kernel
        filter_gaussian_kernel<<<gridSize, blockSize>>>(d_xyz, d_val, d_ar1, d_ar2, d_res, radius, sigma, values.size(), grid_x, grid_y, grid_z);

        // download result
        DynamicVector result(values.size());
        d_res.download(result.data());
        return result;
    }
#endif
    if(device == solver::CPU){
        logging::error(false, "not supported");
    }
}
