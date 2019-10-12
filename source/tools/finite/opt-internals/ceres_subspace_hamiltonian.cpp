//
// Created by david on 2019-10-09.
//

#include <tools/finite/opt.h>
#include <state/class_finite_state.h>
#include <general/nmspc_omp.h>

Eigen::Tensor<std::complex<double>,6> tools::finite::opt::internal::get_multi_hamiltonian(const class_finite_state & state){
//    if(cache.multiham) return cache.multiham.value();
    auto mpo = state.get_multimpo();
    tools::log->trace("Contracting multi hamiltonian...");
    auto & envL = state.get_ENVL(state.active_sites.front());
    auto & envR = state.get_ENVR(state.active_sites.back());
    if (envL.get_position() != state.active_sites.front()) throw std::runtime_error(fmt::format("Mismatch in ENVL and active site positions: {} != {}", envL.get_position() , state.active_sites.front()));
    if (envR.get_position() != state.active_sites.back())  throw std::runtime_error(fmt::format("Mismatch in ENVR and active site positions: {} != {}", envR.get_position() , state.active_sites.back()));
//    cache.multiham =
    long dim0 = mpo.dimension(2);
    long dim1 = envL.block.dimension(0);
    long dim2 = envR.block.dimension(0);
    Eigen::Tensor<std::complex<double>,6> multiham(dim0,dim1,dim2,dim0,dim1,dim2);
    multiham.device(omp::dev) =
            envL.block
                    .contract(mpo           , Textra::idx({2},{0}))
                    .contract(envR.block    , Textra::idx({2},{2}))
                    .shuffle(Textra::array6{2,0,4,3,1,5});
    tools::log->trace("Contracting multi hamiltonian... OK");
//    return cache.multiham.value();
    return multiham;
}

Eigen::Tensor<std::complex<double>,6>   tools::finite::opt::internal::get_multi_hamiltonian2(const class_finite_state & state) {
//    if(cache.multiham_sq) return cache.multiham_sq.value();
    auto mpo = state.get_multimpo();
    tools::log->trace("Contracting multi hamiltonian squared...");
    auto & env2L = state.get_ENV2L(state.active_sites.front());
    auto & env2R = state.get_ENV2R(state.active_sites.back());
    if (env2L.get_position() != state.active_sites.front()) throw std::runtime_error(fmt::format("Mismatch in ENVL and active site positions: {} != {}", env2L.get_position() , state.active_sites.front()));
    if (env2R.get_position() != state.active_sites.back())  throw std::runtime_error(fmt::format("Mismatch in ENVR and active site positions: {} != {}", env2R.get_position() , state.active_sites.back()));
    long dim0 = mpo.dimension(2);
    long dim1 = env2L.block.dimension(0);
    long dim2 = env2R.block.dimension(0);
//    cache.multiham_sq =
    Eigen::Tensor<std::complex<double>,6> multiham_sq(dim0,dim1,dim2,dim0,dim1,dim2);
    multiham_sq.device(omp::dev) =
            env2L.block
                    .contract(mpo             , Textra::idx({2},{0}))
                    .contract(mpo             , Textra::idx({5,2},{2,0}))
                    .contract(env2R.block     , Textra::idx({2,4},{2,3}))
                    .shuffle(Textra::array6{2,0,4,3,1,5});
    tools::log->trace("Contracting multi hamiltonian squared... OK");
    return multiham_sq;
//    return cache.multiham_sq.value();
}

Eigen::MatrixXcd tools::finite::opt::internal::get_multi_hamiltonian_matrix(const class_finite_state & state) {
//    if(cache.multiham_mat) return cache.multiham_mat.value();
    long size = state.active_problem_size();
    auto ham_tensor = tools::finite::opt::internal::get_multi_hamiltonian(state);
    auto cols       = ham_tensor.dimension(0)* ham_tensor.dimension(1)* ham_tensor.dimension(2);
    auto rows       = ham_tensor.dimension(3)* ham_tensor.dimension(4)* ham_tensor.dimension(5);

    if(rows != size)
        throw std::runtime_error (fmt::format("Mismatch hamiltonian dim0*dim1*dim2 and cols: {} != {}",cols, size));
    if(cols != size)
        throw std::runtime_error (fmt::format("Mismatch hamiltonian dim3*dim4*dim5 and rows: {} != {}",rows, size));
    return Eigen::Map<Eigen::MatrixXcd> (ham_tensor.data(),size,size).transpose();
//    cache.multiham_mat =  Eigen::Map<MType> (ham_tensor.data(),size,size).transpose();
//    return cache.multiham_mat.value();
}


Eigen::MatrixXcd tools::finite::opt::internal::get_multi_hamiltonian2_matrix(const class_finite_state & state) {
//    if(cache.multiham_sq_mat) return cache.multiham_sq_mat.value();
    long size = state.active_problem_size();
    auto ham_squared_tensor = tools::finite::opt::internal::get_multi_hamiltonian2(state);
    auto cols       = ham_squared_tensor.dimension(0)* ham_squared_tensor.dimension(1)* ham_squared_tensor.dimension(2);
    auto rows       = ham_squared_tensor.dimension(3)* ham_squared_tensor.dimension(4)* ham_squared_tensor.dimension(5);
    if(rows != size)
        throw std::runtime_error (fmt::format("Mismatch hamiltonian sq dim0*dim1*dim2 and cols: {} != {}",cols, size));
    if(cols != size)
        throw std::runtime_error (fmt::format("Mismatch hamiltonian sq dim3*dim4*dim5 and rows: {} != {}",rows, size));
    return Eigen::Map<Eigen::MatrixXcd> (ham_squared_tensor.data(),size,size).transpose();
//    cache.multiham_sq_mat  = Eigen::Map<MType> (ham_squared_tensor.data(),size,size).transpose();
//    return cache.multiham_sq_mat.value();
}


Eigen::MatrixXcd  tools::finite::opt::internal::get_multi_hamiltonian2_subspace_matrix_new(const class_finite_state & state,const Eigen::MatrixXcd & eigvecs ){
//    if(cache.multiham_sq_sub) return cache.multiham_sq_sub.value();
    auto mpo = state.get_multimpo();
    auto & env2L = state.get_ENV2L(state.active_sites.front()).block;
    auto & env2R = state.get_ENV2R(state.active_sites.back()).block;

    tools::log->trace("Contracting hamiltonian squared matrix in subspace new...");
    auto dims = state.active_dimensions();
    size_t log2chiL  = std::log2(dims[1]);
    size_t log2chiR  = std::log2(dims[2]);
    size_t log2spin  = std::log2(dims[0]);
    size_t eignum    = eigvecs.cols(); //Number of eigenvectors
    size_t eigdim    = eigvecs.rows(); //Length of each eigenvector
    using map        = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>;

    Eigen::Tensor<Scalar,0> H2_ij;
    Eigen::Tensor<Scalar,3> Hv(dims);
    Eigen::MatrixXcd H2(eignum,eignum);


    if(log2spin > log2chiL + log2chiR){
        if (log2chiL >= log2chiR){
            tools::log->trace("get_H2 path: log2spin > log2chiL + log2chiR  and  log2chiL >= log2chiR");
            for (size_t col = 0; col < eignum; col++ ){
                auto theta_j = map(eigvecs.data() + col*eigdim, dims);
                Hv.device(omp::dev) =
                    theta_j
                     .contract(env2L , Textra::idx({1}, {0}))
                     .contract(mpo   , Textra::idx({0,3}, {2,0}))
                     .contract(env2R , Textra::idx({0,3}, {0,2}))
                     .contract(mpo   , Textra::idx({2,1,4}, {2,0,1}))
                     .shuffle(         Textra::array3{2,0,1});
                for (size_t row = col; row < eignum; row++ ){
                    auto theta_i = map(eigvecs.data()+row*eigdim, dims);
                    H2_ij.device(omp::dev) = theta_i.conjugate().contract(Hv, Textra::idx({0,1,2},{0,1,2}));
                    H2(row,col) = H2_ij(0);
                }
            }
        }
        else{
            tools::log->trace("get_H2 path: log2spin > log2chiL + log2chiR  and  log2chiL < log2chiR");
            for (size_t col = 0; col < eignum; col++ ){
                auto theta_j = map(eigvecs.data() + col*eigdim, dims);
                Hv.device(omp::dev) =
                    theta_j
                     .contract(env2R    , Textra::idx({2}, {0}))
                     .contract(mpo      , Textra::idx({0,3}, {2,1}))
                     .contract(env2L    , Textra::idx({0,3}, {0,2}))
                     .contract(mpo      , Textra::idx({2,4,1}, {2,0,1}))
                     .shuffle(            Textra::array3{2,1,0});
                for (size_t row = col; row < eignum; row++ ){
                    auto theta_i = map(eigvecs.data()+ row*eigdim, dims);
                    H2_ij.device(omp::dev) = theta_i.conjugate().contract(Hv, Textra::idx({0,1,2},{0,1,2}));
                    H2(row,col) = H2_ij(0);
                }
            }
        }
    }else{
        tools::log->trace("get_H2 path: log2spin <= log2chiL + log2chiR");
        for (size_t col = 0; col < eignum; col++ ){
            auto theta_j = map(eigvecs.data()+ col*eigdim, dims);
            Hv.device(omp::dev) =
                theta_j
                 .contract(env2L , Textra::idx({1}, {0}))
                 .contract(mpo   , Textra::idx({0,3}, {2,0}))
                 .contract(mpo   , Textra::idx({4,2}, {2,0}))
                 .contract(env2R , Textra::idx({0,2,3}, {0,2,3}))
                 .shuffle(         Textra::array3{1,0,2});
            for (size_t row = col; row < eignum; row++ ){
                auto theta_i = map(eigvecs.data() + row*eigdim, dims);
                H2_ij.device(omp::dev) = theta_i.conjugate().contract(Hv, Textra::idx({0,1,2},{0,1,2}));
                H2(row,col) = H2_ij(0);
            }
        }
    }

    H2 = H2.selfadjointView<Eigen::Lower>();
    tools::log->trace("Contracting hamiltonian squared matrix in subspace new... OK");

    if(H2.hasNaN())throw std::runtime_error("H2 has NaN's!");
//    Eigen::MatrixXcd H2_old = tools::finite::opt::internal::get_multi_hamiltonian2_subspace_matrix(state, eigvecs);
//
//    if(not H2.isApprox(H2_old,1e-4)){
//        std::cout << "H2 new = \n" << H2 << std::endl;
//        std::cout << "H2 old = \n" << H2_old << std::endl;
//        tools::log->warn("H2 subspace mismatch: {:.16f}", (H2 - H2_old).cwiseAbs().sum());
//    }
    if(not H2.isApprox(H2.adjoint(), H2.size()*1e-12)){
        throw std::runtime_error(fmt::format("H2_subspace is not hermitian: {:.16f}", (H2 - H2.adjoint()).cwiseAbs().sum()));
    }
//    return H2.conjugate();
    return H2;

}



Eigen::MatrixXcd  tools::finite::opt::internal::get_multi_hamiltonian2_subspace_matrix(const class_finite_state & state,const Eigen::MatrixXcd & eigvecs ){
//    if(cache.multiham_sq_sub) return cache.multiham_sq_sub.value();
    auto mpo = state.get_multimpo();
    tools::log->trace("Contracting hamiltonian squared matrix in subspace old...");
    auto dims = state.active_dimensions();
    Eigen::DSizes<long,4> eigvecs_dims {dims[0],dims[1],dims[2],eigvecs.cols()};
    auto eigvecs_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(eigvecs.data(), eigvecs_dims );

    auto & env2L = state.get_ENV2L(state.active_sites.front()).block;
    auto & env2R = state.get_ENV2R(state.active_sites.back()).block;

    size_t log2chiL  = std::log2(dims[1]);
    size_t log2chiR  = std::log2(dims[2]);
    size_t log2spin  = std::log2(dims[0]);
    long dimH2 = eigvecs.cols();
    Eigen::Tensor<Scalar,2> H2(dimH2,dimH2);
    if(log2spin > log2chiL + log2chiR){
        if (log2chiL >= log2chiR){
            tools::log->trace("get_H2 path: log2spin > log2chiL + log2chiR  and  log2chiL >= log2chiR ");
            H2 =
                eigvecs_tensor
                        .contract(env2L,                      Textra::idx({1},{0}))
                        .contract(mpo  ,                      Textra::idx({0,4},{2,0}))
                        .contract(env2R,                      Textra::idx({0,4},{0,2}))
                        .contract(mpo  ,                      Textra::idx({3,2,5},{2,0,1}))
                        .contract(eigvecs_tensor.conjugate(), Textra::idx({3,1,2},{0,1,2}))
                        .shuffle(                             Textra::array2{1,0});
        }
        else{
            tools::log->trace("get_H2 path: log2spin > log2chiL + log2chiR  and  log2chiL < log2chiR ");
            H2 =
                eigvecs_tensor
                        .contract(env2R,                      Textra::idx({2},{0}))
                        .contract(mpo  ,                      Textra::idx({0,4},{2,1}))
                        .contract(env2L,                      Textra::idx({0,4},{0,2}))
                        .contract(mpo  ,                      Textra::idx({3,2,5},{2,1,0}))
                        .contract(eigvecs_tensor.conjugate(), Textra::idx({3,2,1},{0,1,2}))
                        .shuffle(                             Textra::array2{1,0});
        }
    }else{
        tools::log->trace("get_H2 path: log2spin <= log2chiL + log2chiR");
        H2 =
//            eigvecs_tensor.conjugate()
//                    .contract(env2L,                      Textra::idx({1},{1}))
//                    .contract(mpo  ,                      Textra::idx({0,5},{3,0}))
//                    .contract(mpo  ,                      Textra::idx({5,3},{3,0}))
//                    .contract(eigvecs_tensor,             Textra::idx({5,2},{0,1}))
//                    .contract(env2R,                      Textra::idx({4,0,3,2},{0,1,2,3}));

            eigvecs_tensor
                    .contract(env2L,                      Textra::idx({1},{0}))
                    .contract(mpo  ,                      Textra::idx({0,4},{2,0}))
                    .contract(mpo  ,                      Textra::idx({5,3},{2,0}))
                    .contract(eigvecs_tensor.conjugate(), Textra::idx({5,2},{0,1}))
                    .contract(env2R,                      Textra::idx({0,4,2,3},{0,1,2,3}))
                    .shuffle(                             Textra::array2{1,0});
    }
    tools::log->trace("Contracting hamiltonian squared matrix in subspace old... OK");
    return Eigen::Map<Eigen::MatrixXcd>(H2.data(),H2.dimension(0),H2.dimension(1));
}



