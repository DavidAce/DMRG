


#include <Eigen/Core>
#include <h5pp/h5pp.h>
#include <io/nmspc_filesystem.h>
#include <math/class_svd_wrapper.h>

int main(){
    using reciter = fs::recursive_directory_iterator;
    for (auto & item : reciter(std::string(TEST_MATRIX_DIR))){
        if(item.path().filename().string().find("svdmatrix") == std::string::npos) continue;
        if(item.path().extension() != ".h5") continue;
        std::cout <<  "item: " << item << std::endl;
        int logLevel = 1;
        h5pp::File file(item.path().string(), h5pp::FilePermission::READONLY, logLevel);
        Eigen::MatrixXcd matrix;
        file.readDataset(matrix,"svdmatrix");
        class_SVD SVD;
        auto [U1,S1,V1] = SVD.decompose(matrix);
        SVD.use_lapacke = true;
        auto [U2,S2,V2] = SVD.decompose(matrix);

        std::cout << "Abs diff S2-S1: " << (S2.array()-S1.array()).cwiseAbs().sum() << std::endl;
    }

    return 0;

}