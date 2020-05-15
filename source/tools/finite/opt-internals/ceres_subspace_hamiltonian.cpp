//
// Created by david on 2019-10-09.
//

#include <tools/finite/opt.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <config/nmspc_settings.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

Eigen::Tensor<std::complex<double>,6> tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_tensor(const class_model_finite & model, const class_edges_finite & edges){
    const auto & mpo = model.get_multisite_mpo();
    const auto & env = edges.get_multisite_ene_blk();
    tools::log->trace("Contracting multisite hamiltonian");
    tools::common::profile::t_ham->tic();

    long dim0 = mpo.dimension(2);
    long dim1 = env.L.dimension(0);
    long dim2 = env.R.dimension(0);
    OMP omp(settings::threading::num_threads);
    Eigen::Tensor<std::complex<double>,6> ham(dim0, dim1, dim2, dim0, dim1, dim2);
    ham.device(omp.dev) =
            env.L
                    .contract(mpo      , Textra::idx({2},{0}))
                    .contract(env.R    , Textra::idx({2},{2}))
                    .shuffle(Textra::array6{2,0,4,3,1,5});
    auto cols       = ham.dimension(0)* ham.dimension(1)* ham.dimension(2);
    auto rows       = ham.dimension(3)* ham.dimension(4)* ham.dimension(5);
    long size       = dim0*dim1*dim2;
    if(rows != size) throw std::runtime_error (fmt::format("Mismatch in multisite hamiltonian dim0*dim1*dim2 and cols: {} != {}",cols, size));
    if(cols != size) throw std::runtime_error (fmt::format("Mismatch in multisite hamiltonian dim3*dim4*dim5 and rows: {} != {}",rows, size));

    auto ham_map = Eigen::Map<Eigen::MatrixXcd>(ham.data(), rows,cols);
    double non_hermiticity = (ham_map - ham_map.adjoint()).cwiseAbs().sum()/static_cast<double>(ham_map.size());
    double sparcity = static_cast<double>((ham_map.array().cwiseAbs2() != 0.0).count())/static_cast<double>(ham_map.size());

    if(non_hermiticity > 1e-12) throw std::runtime_error(fmt::format("multisite hamiltonian is not hermitian: {:.16f}",non_hermiticity));
    if(non_hermiticity > 1e-14) tools::log->warn("multisite hamiltonian is slightly non-hermitian: {:.16f}",non_hermiticity);
    if(ham_map.hasNaN())        throw std::runtime_error("multisite hamiltonian has NaN's!");
    tools::log->trace("multisite hamiltonian nonzeros: {:.8f} %", sparcity*100);
    tools::common::profile::t_ham->toc();
    return ham;
}

Eigen::Tensor<std::complex<double>,6>   tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_tensor(const class_model_finite & model, const class_edges_finite & edges) {
    const auto & mpo = model.get_multisite_mpo();
    const auto & env = edges.get_multisite_var_blk();

    tools::log->trace("Contracting multisite hamiltonian squared");
    tools::common::profile::t_hsq->tic();

    long dim0 = mpo.dimension(2);
    long dim1 = env.L.dimension(0);
    long dim2 = env.R.dimension(0);
    OMP omp(settings::threading::num_threads);
    Eigen::Tensor<std::complex<double>,6> ham_sq(dim0, dim1, dim2, dim0, dim1, dim2);
    ham_sq.device(omp.dev) =
            env.L
                    .contract(mpo             , Textra::idx({2},{0}))
                    .contract(mpo             , Textra::idx({5,2},{2,0}))
                    .contract(env.R     , Textra::idx({2,4},{2,3}))
                    .shuffle(Textra::array6{2,0,4,3,1,5});
    auto cols       = ham_sq.dimension(0)* ham_sq.dimension(1)* ham_sq.dimension(2);
    auto rows       = ham_sq.dimension(3)* ham_sq.dimension(4)* ham_sq.dimension(5);
    long size       = dim0*dim1*dim2;
    if(rows != size) throw std::runtime_error (fmt::format("Mismatch in multisite hamiltonian squared dim0*dim1*dim2 and cols: {} != {}",cols, size));
    if(cols != size) throw std::runtime_error (fmt::format("Mismatch in multisite hamiltonian squared dim3*dim4*dim5 and rows: {} != {}",rows, size));

    auto ham_sq_map = Eigen::Map<Eigen::MatrixXcd>(ham_sq.data(), rows,cols);
    double non_hermiticity = (ham_sq_map - ham_sq_map.adjoint()).cwiseAbs().sum()/static_cast<double>(ham_sq_map.size());
    double sparcity = static_cast<double>((ham_sq_map.array().cwiseAbs2() != 0.0).count())/static_cast<double>(ham_sq_map.size());

    if(non_hermiticity > 1e-12) throw std::runtime_error(fmt::format("multisite hamiltonian squared is not hermitian: {:.16f}",non_hermiticity));
    if(non_hermiticity > 1e-14) tools::log->warn("multisite hamiltonian squared is slightly non-hermitian: {:.16f}",non_hermiticity);
    if(ham_sq_map.hasNaN())     throw std::runtime_error("multisite hamiltonian squared has NaN's!");
    tools::log->trace("multisite hamiltonian squared nonzeros: {:.8f} %", sparcity*100);
    tools::common::profile::t_hsq->toc();
    return ham_sq;
}



Eigen::MatrixXcd tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_matrix(const class_model_finite & model, const class_edges_finite & edges) {
    auto ham_tensor = tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_tensor(model,edges);
    auto rows = ham_tensor.dimension(0);
    auto cols = ham_tensor.dimension(1);
    return Eigen::Map<Eigen::MatrixXcd> (ham_tensor.data(),rows,cols).transpose().selfadjointView<Eigen::Lower>();
}


Eigen::MatrixXcd tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_matrix(const class_model_finite & model, const class_edges_finite & edges) {
    auto ham_squared_tensor = tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_tensor(model,edges);
    auto rows = ham_squared_tensor.dimension(0);
    auto cols = ham_squared_tensor.dimension(1);
    return Eigen::Map<Eigen::MatrixXcd> (ham_squared_tensor.data(),rows,cols).transpose().selfadjointView<Eigen::Lower>();
}


Eigen::MatrixXcd  tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs ){
    const auto & mpo = model.get_multisite_mpo();
    const auto & env = edges.get_multisite_var_blk();
    tools::common::profile::t_hsq->tic();
    tools::log->trace("Contracting subspace hamiltonian squared new");
    long dim0 = mpo.dimension(2);
    long dim1 = env.L.dimension(0);
    long dim2 = env.R.dimension(0);
    double log2spin  = std::log2(dim0);
    double log2chiL  = std::log2(dim1);
    double log2chiR  = std::log2(dim2);
    long   eignum    = eigvecs.cols(); //Number of eigenvectors
    long   eigdim    = eigvecs.rows(); //Length of each eigenvector
    using map        = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>;

    Eigen::Tensor<Scalar,0> H2_ij;
    Eigen::Tensor<Scalar,3> Hv(dim0,dim1,dim2);
    Eigen::MatrixXcd H2(eignum,eignum);

    OMP omp(settings::threading::num_threads);

    if(log2spin >= std::max(log2chiL,log2chiR)){
        if (log2chiL >= log2chiR){
            tools::log->trace("get_H2 path: log2spin >= std::max(log2chiL, log2chiR)  and  log2chiL >= log2chiR");
            for (auto col = 0; col < eignum; col++ ){
                auto theta_j = map(eigvecs.data() + col*eigdim, dim0,dim1,dim2);
                Hv.device(omp.dev) =
                    theta_j
                     .contract(env.L , Textra::idx({1}, {0}))
                     .contract(mpo   , Textra::idx({0,3}, {2,0}))
                     .contract(env.R , Textra::idx({0,3}, {0,2}))
                     .contract(mpo   , Textra::idx({2,1,4}, {2,0,1}))
                     .shuffle(         Textra::array3{2,0,1});
                for (auto row = col; row < eignum; row++ ){
                    auto theta_i = map(eigvecs.data()+row*eigdim, dim0,dim1,dim2);
                    H2_ij.device(omp.dev) = theta_i.conjugate().contract(Hv, Textra::idx({0,1,2},{0,1,2}));
                    H2(row,col) = H2_ij(0);
                }
            }
        }
        else{
            tools::log->trace("get_H2 path: log2spin >= std::max(log2chiL, log2chiR)  and  log2chiL < log2chiR");
            for (auto col = 0; col < eignum; col++ ){
                auto theta_j = map(eigvecs.data() + col*eigdim, dim0,dim1,dim2);
                Hv.device(omp.dev) =
                    theta_j
                     .contract(env.R    , Textra::idx({2}, {0}))
                     .contract(mpo      , Textra::idx({0,3}, {2,1}))
                     .contract(env.L    , Textra::idx({0,3}, {0,2}))
                     .contract(mpo      , Textra::idx({2,4,1}, {2,0,1}))
                     .shuffle(            Textra::array3{2,1,0});
                for (auto row = col; row < eignum; row++ ){
                    auto theta_i = map(eigvecs.data()+ row*eigdim, dim0,dim1,dim2);
                    H2_ij.device(omp.dev) = theta_i.conjugate().contract(Hv, Textra::idx({0,1,2},{0,1,2}));
                    H2(row,col) = H2_ij(0);
                }
            }
        }
    }else{
        tools::log->trace("get_H2 path: log2spin <= log2chiL + log2chiR");
        for (auto col = 0; col < eignum; col++ ){
            auto theta_j = map(eigvecs.data()+ col*eigdim, dim0,dim1,dim2);
            Hv.device(omp.dev) =
                theta_j
                 .contract(env.L , Textra::idx({1}, {0}))
                 .contract(mpo   , Textra::idx({0,3}, {2,0}))
                 .contract(mpo   , Textra::idx({4,2}, {2,0}))
                 .contract(env.R , Textra::idx({0,2,3}, {0,2,3}))
                 .shuffle(         Textra::array3{1,0,2});
            for (auto row = col; row < eignum; row++ ){
                auto theta_i = map(eigvecs.data() + row*eigdim, dim0,dim1,dim2);
                H2_ij.device(omp.dev) = theta_i.conjugate().contract(Hv, Textra::idx({0,1,2},{0,1,2}));
                H2(row,col) = H2_ij(0);
            }
        }
    }

    H2 = H2.selfadjointView<Eigen::Lower>();
    double non_hermiticity = (H2 - H2.adjoint()).cwiseAbs().sum()/static_cast<double>(H2.size());
    double sparcity = static_cast<double>((H2.array().cwiseAbs2() != 0.0).count())/static_cast<double>(H2.size());

    if(non_hermiticity > 1e-12) throw std::runtime_error(fmt::format("subspace hamiltonian squared is not hermitian: {:.16f}",non_hermiticity));
    if(non_hermiticity > 1e-14) tools::log->warn("subspace hamiltonian squared is slightly non-hermitian: {:.16f}",non_hermiticity);
    if(H2.hasNaN())     throw std::runtime_error("subspace hamiltonian squared has NaN's!");
    tools::log->trace("multisite hamiltonian squared nonzeros: {:.8f} %", sparcity*100);
    tools::common::profile::t_hsq->toc();
    return H2;

}



Eigen::MatrixXcd  tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix(const class_model_finite & model,const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs ){
//    if(cache.multiham_sq_sub) return cache.multiham_sq_sub.value();
    const auto & mpo = model.get_multisite_mpo();
    const auto & env = edges.get_multisite_var_blk();
    tools::log->trace("Contracting hamiltonian squared matrix in subspace old");


    long dim0 = mpo.dimension(2);
    long dim1 = env.L.dimension(0);
    long dim2 = env.R.dimension(0);
    double log2spin  = std::log2(dim0);
    double log2chiL  = std::log2(dim1);
    double log2chiR  = std::log2(dim2);
    long dimH2 = eigvecs.cols();
    Eigen::DSizes<long,4> eigvecs_dims {dim0,dim1,dim2,eigvecs.cols()};
    auto eigvecs_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(eigvecs.data(), eigvecs_dims );
    Eigen::Tensor<Scalar,2> H2(dimH2,dimH2);

    if(log2spin >= std::max(log2chiL, log2chiR)){
        if (log2chiL >= log2chiR){
            tools::log->trace("get_H2 path: log2spin >= std::max(log2chiL, log2chiR)  and  log2chiL >= log2chiR ");
            H2 =
                eigvecs_tensor
                        .contract(env.L,                      Textra::idx({1},{0}))
                        .contract(mpo  ,                      Textra::idx({0,4},{2,0}))
                        .contract(env.R,                      Textra::idx({0,4},{0,2}))
                        .contract(mpo  ,                      Textra::idx({3,2,5},{2,0,1}))
                        .contract(eigvecs_tensor.conjugate(), Textra::idx({3,1,2},{0,1,2}))
                        .shuffle(                             Textra::array2{1,0});
        }
        else{
            tools::log->trace("get_H2 path: log2spin >= std::max(log2chiL, log2chiR)  and  log2chiL < log2chiR ");
            H2 =
                eigvecs_tensor
                        .contract(env.R,                      Textra::idx({2},{0}))
                        .contract(mpo  ,                      Textra::idx({0,4},{2,1}))
                        .contract(env.L,                      Textra::idx({0,4},{0,2}))
                        .contract(mpo  ,                      Textra::idx({3,2,5},{2,1,0}))
                        .contract(eigvecs_tensor.conjugate(), Textra::idx({3,2,1},{0,1,2}))
                        .shuffle(                             Textra::array2{1,0});
        }
    }else{
        tools::log->trace("get_H2 path: log2spin < std::max(log2chiL, log2chiR)");
        H2 =
//            eigvecs_tensor.conjugate()
//                    .contract(env2L,                      Textra::idx({1},{1}))
//                    .contract(mpo  ,                      Textra::idx({0,5},{3,0}))
//                    .contract(mpo  ,                      Textra::idx({5,3},{3,0}))
//                    .contract(eigvecs_tensor,             Textra::idx({5,2},{0,1}))
//                    .contract(env2R,                      Textra::idx({4,0,3,2},{0,1,2,3}));

            eigvecs_tensor
                    .contract(env.L,                      Textra::idx({1},{0}))
                    .contract(mpo  ,                      Textra::idx({0,4},{2,0}))
                    .contract(mpo  ,                      Textra::idx({5,3},{2,0}))
                    .contract(eigvecs_tensor.conjugate(), Textra::idx({5,2},{0,1}))
                    .contract(env.R,                      Textra::idx({0,4,2,3},{0,1,2,3}))
                    .shuffle(                             Textra::array2{1,0});
    }
    auto H2_map = Eigen::Map<Eigen::MatrixXcd>(H2.data(),H2.dimension(0),H2.dimension(1));
    double non_hermiticity = (H2_map - H2_map.adjoint()).cwiseAbs().sum()/static_cast<double>(H2.size());
    double sparcity = static_cast<double>((H2_map.array().cwiseAbs2() != 0.0).count())/static_cast<double>(H2.size());

    if(non_hermiticity > 1e-12) throw std::runtime_error(fmt::format("subspace hamiltonian squared is not hermitian: {:.16f}",non_hermiticity));
    if(non_hermiticity > 1e-14) tools::log->warn("subspace hamiltonian squared is slightly non-hermitian: {:.16f}",non_hermiticity);
    if(H2_map.hasNaN())         throw std::runtime_error("subspace hamiltonian squared has NaN's!");
    tools::log->trace("multisite hamiltonian squared nonzeros: {:.8f} %", sparcity*100);


    return H2_map;
}







// Template definitions

using Scalar = std::complex<double>;
template<typename T> using MatrixType = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;


template <typename T>
Eigen::Tensor<T,6> tools::finite::opt::internal::get_multi_hamiltonian_tensor(const class_model_finite & model, const class_edges_finite & edges){
    static_assert(std::is_same<T,std::complex<double>>::value or std::is_same<T,double>::value,"Wrong type, expected double or complex double");
    if      constexpr(std::is_same<T,std::complex<double>>::value) return tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_tensor(model,edges);
    else if constexpr(std::is_same<T,double>::value)               return tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_tensor(model,edges).real();
}

// Explicit instantiations
template Eigen::Tensor<double,6>tools::finite::opt::internal::get_multi_hamiltonian_tensor<double>(const class_model_finite & model, const class_edges_finite & edges);
template Eigen::Tensor<Scalar,6>tools::finite::opt::internal::get_multi_hamiltonian_tensor<Scalar>(const class_model_finite & model, const class_edges_finite & edges);


template <typename T>
Eigen::Tensor<T,6> tools::finite::opt::internal::get_multi_hamiltonian_squared_tensor(const class_model_finite & model, const class_edges_finite & edges){
    static_assert(std::is_same<T,std::complex<double>>::value or std::is_same<T,double>::value,"Wrong type, expected double or complex double");
    if      constexpr(std::is_same<T,std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_squared_tensor(model,edges);
    else if constexpr(std::is_same<T,double>::value)               return local_hamiltonians::get_multi_hamiltonian_squared_tensor(model,edges).real();
}
// Explicit instantiations
template Eigen::Tensor<double,6> tools::finite::opt::internal::get_multi_hamiltonian_squared_tensor<double>(const class_model_finite & model, const class_edges_finite & edges);
template Eigen::Tensor<Scalar,6> tools::finite::opt::internal::get_multi_hamiltonian_squared_tensor<Scalar>(const class_model_finite & model, const class_edges_finite & edges);


template <typename T>
MatrixType<T> tools::finite::opt::internal::get_multi_hamiltonian_matrix(const class_model_finite & model, const class_edges_finite & edges){
    static_assert(std::is_same<T,std::complex<double>>::value or std::is_same<T,double>::value,"Wrong type, expected double or complex double");
    if      constexpr(std::is_same<T,std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_matrix(model,edges);
    else if constexpr(std::is_same<T,double>::value)               return local_hamiltonians::get_multi_hamiltonian_matrix(model,edges).real();
}
// Explicit instantiations
template MatrixType<double> tools::finite::opt::internal::get_multi_hamiltonian_matrix<double>(const class_model_finite & model, const class_edges_finite & edges);
template MatrixType<Scalar> tools::finite::opt::internal::get_multi_hamiltonian_matrix<Scalar>(const class_model_finite & model, const class_edges_finite & edges);

template <typename T>
MatrixType<T> tools::finite::opt::internal::get_multi_hamiltonian_squared_matrix(const class_model_finite & model, const class_edges_finite & edges){
    static_assert(std::is_same<T,std::complex<double>>::value or std::is_same<T,double>::value,"Wrong type, expected double or complex double");
    if      constexpr(std::is_same<T,std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_squared_matrix(model,edges);
    else if constexpr(std::is_same<T,double>::value)               return local_hamiltonians::get_multi_hamiltonian_squared_matrix(model,edges).real();
}
// Explicit instantiations
template MatrixType<double> tools::finite::opt::internal::get_multi_hamiltonian_squared_matrix<double>(const class_model_finite & model, const class_edges_finite & edges);
template MatrixType<Scalar> tools::finite::opt::internal::get_multi_hamiltonian_squared_matrix<Scalar>(const class_model_finite & model, const class_edges_finite & edges);

template <typename T>
MatrixType<T> tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs){
    static_assert(std::is_same<T,std::complex<double>>::value or std::is_same<T,double>::value,"Wrong type, expected double or complex double");
    if      constexpr(std::is_same<T,std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix(model,edges, eigvecs);
    else if constexpr(std::is_same<T,double>::value)               return local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix(model,edges, eigvecs).real();
}
// Explicit instantiations
template MatrixType<double> tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix<double>(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs);
template MatrixType<Scalar> tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix<Scalar>(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs);

template <typename T>
MatrixType<T> tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs){
    static_assert(std::is_same<T,std::complex<double>>::value or std::is_same<T,double>::value,"Wrong type, expected double or complex double");
    if      constexpr(std::is_same<T,std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix_new(model,edges, eigvecs);
    else if constexpr(std::is_same<T,double>::value)               return local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix_new(model,edges, eigvecs).real();
}
// Explicit instantiations
template MatrixType<double> tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new<double>(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs);
template MatrixType<Scalar> tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new<Scalar>(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs);



