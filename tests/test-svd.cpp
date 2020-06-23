#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <Eigen/Core>
#include <h5pp/h5pp.h>
#include <math/svd.h>

TEST_CASE("Singular value decomposition in Eigen and Lapacke", "[svd]") {
    SECTION("Test all svd decompositions") {
        using reciter = h5pp::fs::recursive_directory_iterator;
        svd::solver svd;
        for(auto &item : reciter(std::string(TEST_MATRIX_DIR))) {
            if(item.path().filename().string().find("svdmatrix") == std::string::npos) continue;
            if(item.path().extension() != ".h5") continue;
            size_t     logLevel = 2;
            h5pp::File file(item.path().string(), h5pp::FilePermission::READONLY, logLevel);
            auto       matrix = file.readDataset<Eigen::MatrixXcd>("svdmatrix");
            svd.use_lapacke   = false;
            auto [U1, S1, V1] = svd.decompose(matrix);
            svd.use_lapacke   = true;
            auto [U2, S2, V2] = svd.decompose(matrix);
            double difference = std::log10((S2.array() - S1.array()).cwiseAbs().sum());
            fmt::print("{:<32} diff {:.24f}\n", item.path().filename().string(), difference);
            REQUIRE(difference < -12);
        }
    }
}