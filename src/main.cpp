
#include "core/core.h"
#include "material/material.h"
#include "reader/file.h"
#include "reader/line.h"
#include "reader/reader.h"

#include <Eigen/Sparse>

int main(int argc, char* argv[]) {

    fem::reader::Reader reader{std::string(argv[1])};
    reader.read();

    return 0;
}

