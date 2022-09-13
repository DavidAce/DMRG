#include <string_view>
#include <unsupported/Eigen/CXX11/Tensor>

void print_tensor(const Eigen::Tensor<std::complex<double>,1> & L, std::string_view msg){
    std::printf("%s\n", msg.data());
    for(long i = 0; i < L.size(); i++) std::printf("(%.16f, %.16f)\n",L[i].real(), L[i].imag());
}
void print_tensor(const Eigen::Tensor<std::complex<double>,2> & L, std::string_view msg){
    std::printf("%s\n", msg.data());
    for(long j = 0; j < L.dimension(1); j++)
        for(long i = 0; i < L.dimension(0); i++) std::printf("(%.16f, %.16f)\n",L(i,j).real(), L(i,j).imag());
}

Eigen::Tensor<std::complex<double>,1> broadcast(Eigen::Tensor<std::complex<double>,1> & tensor, std::array<long,1> bcast){
    // Use this function to avoid a bug in Eigen when broadcasting complex tensors of rank 1, with compiler option -mfma
    // See more here https://gitlab.com/libeigen/eigen/-/issues/2351
    std::array<long,2> bcast2 = {bcast[0], 1};
    std::array<long,2> shape2 = {tensor.size(), 1};
    std::array<long,1> shape1 = {tensor.size()*bcast[0]};
    return tensor.reshape(shape2).broadcast(bcast2).reshape(shape1);
}


int main() {
    Eigen::Tensor<std::complex<double>,1> L(1);
    L.setConstant(1.0);

//    print_tensor(L, "L");


    std::array<long,1> bcast = {4};
    Eigen::Tensor<std::complex<double>,1> Lb = broadcast(L, bcast); // Error happens here


    // The line below should print
    //
    //   (1.0000000000000000, 0.0000000000000000)
    //   (1.0000000000000000, 0.0000000000000000)
    //   (1.0000000000000000, 0.0000000000000000)
    //   (1.0000000000000000, 0.0000000000000000)
    //

    print_tensor(Lb, "L.broadcast({4})");

    if(Lb(1) != 1.0) std::printf("Wrong result");

}
