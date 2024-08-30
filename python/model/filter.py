import numpy as np

#import pycuda if possible
try:
    import pycuda.autoinit
    import pycuda.driver as drv
    from pycuda.compiler import SourceModule
    pycuda_available = True
except:
    pycuda_available = False

import numpy as np
import matplotlib.pyplot as plt
import time
from .symmetry import apply_symmetry
from enum import Enum

class FilterFunction(Enum):
    CONSTANT = 1
    GAUSSIAN = 2
    LINEAR   = 3

if pycuda_available:
    mod = SourceModule("""
    __global__ void filter_gaussian_kernel(const float* xyz,
                                      const float* val,
                                      const int* ar1,
                                      const int* ar2,
                                            float* res,
                                      const float rad,
                                      const float sig,
                                      const int n_nodes,
                                      const int blocks_x,
                                      const int blocks_y,
                                      const int blocks_z,
                                      const int filter_func,
                                      const int min_values_to_filter) {
        int n_id = blockIdx.x * blockDim.x + threadIdx.x;
    
        if(n_id >= n_nodes) return;
    
        // extract location
        float x = xyz[3 * n_id + 0];
        float y = xyz[3 * n_id + 1];
        float z = xyz[3 * n_id + 2];
    
        // compute the block location
        int bucket_x = x / rad;
        int bucket_y = y / rad;
        int bucket_z = z / rad;
    
        float weight = 0;
        float value = 0;
        int count = 0;
        for(int bx = max(bucket_x - 1, 0); bx < min(bucket_x + 2, blocks_x); bx++){
            for(int by = max(bucket_y - 1, 0); by < min(bucket_y + 2, blocks_y); by++){
                for(int bz = max(bucket_z - 1, 0); bz < min(bucket_z + 2, blocks_z); bz++){
                    int id_start = ar1[bx * blocks_y * blocks_z + by * blocks_z + bz];
                    int id_end   = ar1[bx * blocks_y * blocks_z + by * blocks_z + bz + 1];
                    for(int i = id_start; i < id_end; i++){
                        int n_id_2 = ar2[i];
                        
                        float x_2 = xyz[3 * n_id_2 + 0];
                        float y_2 = xyz[3 * n_id_2 + 1];
                        float z_2 = xyz[3 * n_id_2 + 2];
    
                        float dist_squared = (x-x_2) * (x-x_2) + (y-y_2) * (y-y_2) + (z-z_2) * (z-z_2);
    
                        if(dist_squared < rad * rad){
                        
                        
                            float wgt;
                            if(filter_func == 1){
                                wgt = 1;
                            }else if(filter_func == 2){
                                wgt = max(0.0f, sqrt(dist_squared) / sig);
                            } else if (filter_func == 3){
                                wgt = exp(- dist_squared / (2 * sig * sig));
                            }
                            count ++;
                            weight += wgt;
                            value  += wgt * val[n_id_2];
                        }
                    }
                }
            }
        }
        if (count < min_values_to_filter){
            weight += (min_values_to_filter - count);
        }
        res[n_id] = value / weight;
    }
    """)


    min_dist_mod = SourceModule("""
    __global__ void min_dist_kernel(
                                      const float* xyz,
                                      const int* ar1,
                                      const int* ar2,
                                            float* res,
                                            float rad,
                                      const int n_nodes,
                                      const int blocks_x,
                                      const int blocks_y,
                                      const int blocks_z) {
        int n_id = blockIdx.x * blockDim.x + threadIdx.x;
    
        if(n_id >= n_nodes) return;
    
        // extract location
        float x = xyz[3 * n_id + 0];
        float y = xyz[3 * n_id + 1];
        float z = xyz[3 * n_id + 2];
    
        // compute the block location
        int bucket_x = x / rad;
        int bucket_y = y / rad;
        int bucket_z = z / rad;
    
        float min_dist = 1e30;
        for(int bx = max(bucket_x - 1, 0); bx < min(bucket_x + 2, blocks_x); bx++){
            for(int by = max(bucket_y - 1, 0); by < min(bucket_y + 2, blocks_y); by++){
                for(int bz = max(bucket_z - 1, 0); bz < min(bucket_z + 2, blocks_z); bz++){
                    int id_start = ar1[bx * blocks_y * blocks_z + by * blocks_z + bz];
                    int id_end   = ar1[bx * blocks_y * blocks_z + by * blocks_z + bz + 1];
                    for(int i = id_start; i < id_end; i++){
                        int n_id_2 = ar2[i];
                        
                        float x_2 = xyz[3 * n_id_2 + 0];
                        float y_2 = xyz[3 * n_id_2 + 1];
                        float z_2 = xyz[3 * n_id_2 + 2];
    
                        float dist_squared = (x-x_2) * (x-x_2) + (y-y_2) * (y-y_2) + (z-z_2) * (z-z_2);
                        if (dist_squared > 1e-24)
                            min_dist = min(dist_squared, min_dist);
                    }
                }
            }
        }
        res[n_id] = sqrt(min_dist);
    }
    """)

class Filter:
    def __init__(self, coords, sigma, symmetries={}, filter_func=FilterFunction.GAUSSIAN):
        self.sigma = sigma
        self.n_values = len(coords)
        self.filter_func = filter_func
        self.symmetries = symmetries
        if pycuda_available and self.sigma > 0.0:
            self.kernel = mod.get_function("filter_gaussian_kernel")
            self.kernel_min_dist = min_dist_mod.get_function("min_dist_kernel")

            # precompute whatever possible
            self.coords        = self._apply_symmetries(coords=coords)
            self.radius        = self._compute_radius()
            self.bounds        = self._compute_range()
            self.grid          = self._compute_grid()
            self.ar1, self.ar2 = self._compute_buckets()
        else:
            pass

    def _compute_radius(self):
        radius = self.sigma
        if self.filter_func == FilterFunction.LINEAR:
            radius = self.sigma
        if self.filter_func == FilterFunction.CONSTANT:
            radius = self.sigma
        if self.filter_func == FilterFunction.GAUSSIAN:
            radius = 3 * self.sigma
        return radius

    def _apply_symmetries(self, coords):
        new_coords, _ = apply_symmetry(coords=coords,
                                                values=np.zeros(self.n_values),
                                                symmetries=self.symmetries)
        return new_coords

    def _compute_range(self):
        min_x, min_y, min_z = np.min(self.coords, axis=0)
        max_x, max_y, max_z = np.max(self.coords, axis=0)

        # step 3, adjust xyz so that min_x and so on are at 0,0,0
        self.coords -= [min_x, min_y, min_z]
        max_x -= min_x
        max_y -= min_y
        max_z -= min_z
        return (max_x, max_y, max_z)

    def _compute_bucket_count(self):
        grid_x = int(np.floor(self.bounds[0] / self.radius)) + 1
        grid_y = int(np.floor(self.bounds[1] / self.radius)) + 1
        grid_z = int(np.floor(self.bounds[2] / self.radius)) + 1
        return (grid_x, grid_y, grid_z)

    def _compute_grid(self):
        while np.prod(self._compute_bucket_count()) > len(self.coords) * 10 or np.prod(self._compute_bucket_count()) < 0:
            self.radius *= 2
        return self._compute_bucket_count()

    def _compute_buckets(self):
        buckets = np.prod(self.grid)

        # step 5, allocate arrays for the grid
        array_1 = np.zeros(buckets + 1, dtype=int)
        array_2 = np.zeros(len(self.coords), dtype=int)

        # step 6, compute how many times a bucket is used
        bucket_node_counts = np.zeros(buckets, dtype=int)
        for m in range(len(self.coords)):
            bucket_x = int(self.coords[m, 0] // self.radius)
            bucket_y = int(self.coords[m, 1] // self.radius)
            bucket_z = int(self.coords[m, 2] // self.radius)
            bucket_idx = bucket_x * self.grid[1] * self.grid[2] \
                       + bucket_y * self.grid[2] \
                       + bucket_z
            bucket_node_counts[bucket_idx] += 1

        # step 7, create a temporary large matrix which holds the indices
        max_nodes_per_bucket = np.max(bucket_node_counts)
        bucket_temp = -np.ones((buckets, max_nodes_per_bucket))
        for m in range(len(self.coords)):
            bucket_x = int(self.coords[m, 0] // self.radius)
            bucket_y = int(self.coords[m, 1] // self.radius)
            bucket_z = int(self.coords[m, 2] // self.radius)
            bucket_idx = bucket_x * self.grid[1] * self.grid[2] \
                       + bucket_y * self.grid[2] \
                       + bucket_z
            bucket_node_counts[bucket_idx] -= 1
            bucket_temp[bucket_idx, bucket_node_counts[bucket_idx]] = m

        # step 8, convert to the previously declared arrays
        for m in range(buckets):
            c = 0
            while c < max_nodes_per_bucket and bucket_temp[m, c] >= 0:
                array_2[array_1[m] + c] = bucket_temp[m, c]
                c += 1
            array_1[m + 1] = array_1[m] + c

        return array_1, array_2

    def apply(self, values):

        if not pycuda_available or self.sigma == 0.0:
            return values

        ghost_factor = round(len(self.coords) / self.n_values)
        values = np.tile(values, ghost_factor)

        # create numpy arrays
        d_xyz = np.array(self.coords  , dtype=np.float32)
        d_val = np.array(values       , dtype=np.float32)
        d_ar1 = np.array(self.ar1     , dtype=np.int32)
        d_ar2 = np.array(self.ar2     , dtype=np.int32)
        d_res = np.zeros(self.n_values, dtype=np.float32)

        block_size = (256, 1, 1)
        grid_size = (int(np.ceil(self.n_values / block_size[0])), 1, 1)

        self.kernel(drv.In(d_xyz),
                    drv.In(d_val),
                    drv.In(d_ar1),
                    drv.In(d_ar2),
                    drv.Out(d_res),
                    np.float32(self.radius ), np.float32(self.sigma  ), np.int32(self.n_values),
                    np.int32  (self.grid[0]), np.int32  (self.grid[1]), np.int32(self.grid[2]  ),
                    np.int32(self.filter_func.value), np.int32(0),
                    block=block_size, grid=grid_size)
        return d_res

    def minimal_distance(self):
        if not pycuda_available or self.sigma == 0.0:
            return 0

        d_xyz = np.array(self.coords  , dtype=np.float32)
        d_ar1 = np.array(self.ar1     , dtype=np.int32)
        d_ar2 = np.array(self.ar2     , dtype=np.int32)
        d_res = np.zeros(self.n_values, dtype=np.float32)

        block_size = (256, 1, 1)
        grid_size  = (int(np.ceil(self.n_values / block_size[0])), 1, 1)

        self.kernel_min_dist(drv.In(d_xyz),
                    drv.In(d_ar1),
                    drv.In(d_ar2),
                    drv.Out(d_res),
                    np.float32(self.radius),
                    np.int32(self.n_values),
                    np.int32(self.grid[0]), np.int32(self.grid[1]), np.int32(self.grid[2]),
                    block=block_size, grid=grid_size)

        return d_res;