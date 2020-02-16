


#include <Eigen/Core>
#include <h5pp/h5pp.h>
#include <io/nmspc_filesystem.h>
#include <math/class_svd_wrapper.h>

int main(){
    namespace fs = tools::fs;
    using reciter = fs::recursive_directory_iterator;
    for (auto & item : reciter(std::string(TEST_DIR) + "/testmatrices")){
        if(item.path().extension() != ".h5") continue;
        std::cout <<  "item: " << item << std::endl;
        int logLevel = 0;
        h5pp::File file(item.path().string(), h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN, logLevel);
        Eigen::MatrixXcd matrix;
        file.readDataset(matrix,"svdmatrix");
        class_SVD SVD;
        auto [U1,S1,V1] = SVD.decompose(matrix);
        SVD.use_lapacke = true;
        auto [U2,S2,V2] = SVD.decompose(matrix);

        std::cout << "Abs diff S2-S1: " << (S2.array()-S1.array()).cwiseAbs().sum() << std::endl;
//        h5pp::File file2(std::string(TEST_DIR) + "/testmatrices_new2/"+ item.path().filename().string(), h5pp::AccessMode::READWRITE, h5pp::CreateMode::TRUNCATE, logLevel);
//        file2.writeDataset(matrix,"svdmatrix",H5D_CONTIGUOUS);
    }

    return 0;

}