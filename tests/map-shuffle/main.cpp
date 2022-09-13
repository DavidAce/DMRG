#include <io/fmt.h>
#include <math/rnd.h>
#include <math/svd.h>
template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 4>, long,long,long,long>
gen_random_tensor(){
    long dL = 2;
    long dR = 2;
    long chiL = rnd::uniform_integer_box(1,8);
    long chiR = rnd::uniform_integer_box(1,8);
    Eigen::Tensor<Scalar, 4> tensor(dL,chiL,dR,chiR);
    tensor.setRandom();
    return std::make_tuple(tensor,dL,chiL,dR,chiR);
}

int main(){
    tools::Logger::setLogger(tools::log,"test",0);


    svd::solver svd;
    for(auto i = 0; i < 100; i++){
        auto [tensor,dL,chiL,dR,chiR] = gen_random_tensor<std::complex<double>>();
        auto [U,S,V] = svd.schmidt(tensor,dL,chiL,dR,chiR);
        if(U.dimension(0) != dL)   throw std::runtime_error(fmt::format("U dim 0 error: {} != {}", U.dimension(0),dL));
        if(U.dimension(1) != chiL) throw std::runtime_error(fmt::format("U dim 1 error: {} != {}", U.dimension(1),chiL));
        if(V.dimension(0) != dR)   throw std::runtime_error(fmt::format("V dim 0 error: {} != {}", V.dimension(0),dR));
        if(V.dimension(2) != chiR) throw std::runtime_error(fmt::format("V dim 2 error: {} != {}", V.dimension(2),chiR));
    }
}