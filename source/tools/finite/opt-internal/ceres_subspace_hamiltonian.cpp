//
// Created by david on 2019-10-09.
//

#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include <config/nmspc_settings.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt.h>
#include <tools/finite/opt_state.h>
using Scalar = std::complex<double>;
using real   = double;
using cplx   = std::complex<double>;
template<typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
MatrixType<T> tools::finite::opt::internal::get_multisite_hamiltonian_matrix(const class_model_finite &model, const class_edges_finite &edges) {
    const auto &mpo = model.get_multisite_tensor();
    const auto &env = edges.get_multisite_ene_blk();
    tools::log->trace("Contracting multisite hamiltonian");
    tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_ham"]->tic();

    long dim0 = mpo.dimension(2);
    long dim1 = env.L.dimension(0);
    long dim2 = env.R.dimension(0);

    Eigen::Tensor<std::complex<double>, 6> ham(dim0, dim1, dim2, dim0, dim1, dim2);
    ham.device(*Textra::omp::dev) = env.L.contract(mpo, Textra::idx({2}, {0})).contract(env.R, Textra::idx({2}, {2})).shuffle(Textra::array6{2, 0, 4, 3, 1, 5});

    auto cols = ham.dimension(0) * ham.dimension(1) * ham.dimension(2);
    auto rows = ham.dimension(3) * ham.dimension(4) * ham.dimension(5);
    long size = dim0 * dim1 * dim2;
    if(rows != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian rows (dim0*dim1*dim2) and size: {} != {}", cols, size));
    if(cols != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian cols (dim3*dim4*dim5) and size: {} != {}", rows, size));

    auto   ham_map         = Eigen::Map<Eigen::MatrixXcd>(ham.data(), rows, cols);
    double non_hermiticity = (ham_map - ham_map.adjoint()).cwiseAbs().sum() / static_cast<double>(ham_map.size());
    double sparcity        = static_cast<double>((ham_map.array().cwiseAbs2() != 0.0).count()) / static_cast<double>(ham_map.size());

    if(non_hermiticity > 1e-12) throw std::runtime_error(fmt::format("multisite hamiltonian is not hermitian: {:.16f}", non_hermiticity));
    if(non_hermiticity > 1e-14) tools::log->warn("multisite hamiltonian is slightly non-hermitian: {:.16f}", non_hermiticity);
    if(ham_map.hasNaN()) throw std::runtime_error("multisite hamiltonian has NaN's!");
    tools::log->trace("multisite hamiltonian nonzeros: {:.8f} %", sparcity * 100);
    tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_ham"]->toc();
    if constexpr(std::is_same_v<T,double>){
        return Eigen::Map<Eigen::MatrixXcd>(ham.data(), rows, cols).real().transpose().selfadjointView<Eigen::Lower>();
    }else{
        return Eigen::Map<Eigen::MatrixXcd>(ham.data(), rows, cols).transpose().selfadjointView<Eigen::Lower>();
    }
}
template MatrixType<real> tools::finite::opt::internal::get_multisite_hamiltonian_matrix(const class_model_finite &model,
                                                                                           const class_edges_finite &edges);
template MatrixType<cplx> tools::finite::opt::internal::get_multisite_hamiltonian_matrix(const class_model_finite &model,
                                                                                           const class_edges_finite &edges);

template<typename T>
MatrixType<T> tools::finite::opt::internal::get_multisite_hamiltonian_squared_subspace_matrix(const class_model_finite &model, const class_edges_finite &edges,
                                                                                              const std::vector<opt_state> &candidate_list) {
    // First, make sure every candidate is actually a basis vector, otherwise this computation would turn difficult if we have to skip rows and columns
    for(const auto &candidate : candidate_list)
        if(not candidate.is_basis_vector)
            throw std::runtime_error("One candidate is not a basis vector. When constructing a hamiltonian subspace matrix, make sure the candidates are all "
                                     "eigenvectors/basis vectors");

    const auto &mpo = model.get_multisite_tensor();
    const auto &env = edges.get_multisite_var_blk();
    tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_hsq"]->tic();
    tools::log->trace("Contracting subspace hamiltonian squared new");
    long   dim0     = mpo.dimension(2);
    long   dim1     = env.L.dimension(0);
    long   dim2     = env.R.dimension(0);
    double log2spin = std::log2(dim0);
    double log2chiL = std::log2(dim1);
    double log2chiR = std::log2(dim2);
    long   eignum   = static_cast<long>(candidate_list.size()); // Number of eigenvectors

    Eigen::Tensor<Scalar, 0> H2_ij;
    Eigen::Tensor<Scalar, 3> Hv(dim0, dim1, dim2);
    Eigen::MatrixXcd         H2(eignum, eignum);

    if(log2spin >= std::max(log2chiL, log2chiR)) {
        if(log2chiL >= log2chiR) {
            tools::log->trace("get_H2 path: log2spin >= std::max(log2chiL, log2chiR)  and  log2chiL >= log2chiR");
            for(auto col = 0; col < eignum; col++) {
                const auto &theta_j = std::next(candidate_list.begin(), col)->get_tensor();
                Hv.device(*Textra::omp::dev)  = theta_j.contract(env.L, Textra::idx({1}, {0}))
                                         .contract(mpo, Textra::idx({0, 3}, {2, 0}))
                                         .contract(env.R, Textra::idx({0, 3}, {0, 2}))
                                         .contract(mpo, Textra::idx({2, 1, 4}, {2, 0, 1}))
                                         .shuffle(Textra::array3{2, 0, 1});
                for(auto row = col; row < eignum; row++) {
                    const auto &theta_i   = std::next(candidate_list.begin(), row)->get_tensor();
                    H2_ij.device(*Textra::omp::dev) = theta_i.conjugate().contract(Hv, Textra::idx({0, 1, 2}, {0, 1, 2}));
                    H2(row, col)          = H2_ij(0);
                }
            }
        } else {
            tools::log->trace("get_H2 path: log2spin >= std::max(log2chiL, log2chiR)  and  log2chiL < log2chiR");
            for(auto col = 0; col < eignum; col++) {
                const auto &theta_j = std::next(candidate_list.begin(), col)->get_tensor();
                Hv.device(*Textra::omp::dev)  = theta_j.contract(env.R, Textra::idx({2}, {0}))
                                         .contract(mpo, Textra::idx({0, 3}, {2, 1}))
                                         .contract(env.L, Textra::idx({0, 3}, {0, 2}))
                                         .contract(mpo, Textra::idx({2, 4, 1}, {2, 0, 1}))
                                         .shuffle(Textra::array3{2, 1, 0});
                for(auto row = col; row < eignum; row++) {
                    const auto &theta_i   = std::next(candidate_list.begin(), row)->get_tensor();
                    H2_ij.device(*Textra::omp::dev) = theta_i.conjugate().contract(Hv, Textra::idx({0, 1, 2}, {0, 1, 2}));
                    H2(row, col)          = H2_ij(0);
                }
            }
        }
    } else {
        tools::log->trace("get_H2 path: log2spin <= log2chiL + log2chiR");
        for(auto col = 0; col < eignum; col++) {
            const auto &theta_j = std::next(candidate_list.begin(), col)->get_tensor();
            Hv.device(*Textra::omp::dev)  = theta_j.contract(env.L, Textra::idx({1}, {0}))
                                     .contract(mpo, Textra::idx({0, 3}, {2, 0}))
                                     .contract(mpo, Textra::idx({4, 2}, {2, 0}))
                                     .contract(env.R, Textra::idx({0, 2, 3}, {0, 2, 3}))
                                     .shuffle(Textra::array3{1, 0, 2});
            for(auto row = col; row < eignum; row++) {
                const auto &theta_i   = std::next(candidate_list.begin(), row)->get_tensor();
                H2_ij.device(*Textra::omp::dev) = theta_i.conjugate().contract(Hv, Textra::idx({0, 1, 2}, {0, 1, 2}));
                H2(row, col)          = H2_ij(0);
            }
        }
    }

    H2                     = H2.selfadjointView<Eigen::Lower>();
    double non_hermiticity = (H2 - H2.adjoint()).cwiseAbs().sum() / static_cast<double>(H2.size());
    double sparcity        = static_cast<double>((H2.array().cwiseAbs2() != 0.0).count()) / static_cast<double>(H2.size());
    if(H2.rows() != H2.cols()) throw std::runtime_error(fmt::format("Hamiltonian tensor is not square: rows [{}] | cols [{}]", H2.rows(), H2.cols()));
    if(non_hermiticity > 1e-12) throw std::runtime_error(fmt::format("subspace hamiltonian squared is not hermitian: {:.16f}", non_hermiticity));
    if(non_hermiticity > 1e-14) tools::log->warn("subspace hamiltonian squared is slightly non-hermitian: {:.16f}", non_hermiticity);
    if(H2.hasNaN()) throw std::runtime_error("subspace hamiltonian squared has NaN's!");
    tools::log->trace("multisite hamiltonian squared nonzeros: {:.8f} %", sparcity * 100);
    tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_hsq"]->toc();
    if constexpr(std::is_same_v<T, double>) return H2.real();
    else
        return H2;
}

// Explicit instantiations
template MatrixType<real> tools::finite::opt::internal::get_multisite_hamiltonian_squared_subspace_matrix(const class_model_finite &           model,
                                                                                                          const class_edges_finite &           edges,
                                                                                                          const std::vector<opt_state> &candidate_list);
template MatrixType<cplx> tools::finite::opt::internal::get_multisite_hamiltonian_squared_subspace_matrix(const class_model_finite &           model,
                                                                                                          const class_edges_finite &           edges,
                                                                                                          const std::vector<opt_state> &candidate_list);

// Eigen::Tensor<std::complex<double>, 6> tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_tensor(const class_model_finite
// &model,
//                                                                                                                              const class_edges_finite &edges)
//                                                                                                                              {
//    const auto &mpo = model.get_multisite_tensor();
//    const auto &env = edges.get_multisite_var_blk();
//
//    tools::log->trace("Contracting multisite hamiltonian squared");
//    tools::common::profile::t_hsq->tic();
//
//    long                                   dim0 = mpo.dimension(2);
//    long                                   dim1 = env.L.dimension(0);
//    long                                   dim2 = env.R.dimension(0);
//    Eigen::Tensor<std::complex<double>, 6> ham_sq(dim0, dim1, dim2, dim0, dim1, dim2);
//    ham_sq.device(*Textra::omp::dev) = env.L.contract(mpo, Textra::idx({2}, {0}))
//                                 .contract(mpo, Textra::idx({5, 2}, {2, 0}))
//                                 .contract(env.R, Textra::idx({2, 4}, {2, 3}))
//                                 .shuffle(Textra::array6{2, 0, 4, 3, 1, 5});
//    auto cols = ham_sq.dimension(0) * ham_sq.dimension(1) * ham_sq.dimension(2);
//    auto rows = ham_sq.dimension(3) * ham_sq.dimension(4) * ham_sq.dimension(5);
//    long size = dim0 * dim1 * dim2;
//    if(rows != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian squared rows (dim0*dim1*dim2) and size: {} != {}", cols, size));
//    if(cols != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian squared cols (dim3*dim4*dim5) and size: {} != {}", rows, size));
//
//    auto   ham_sq_map      = Eigen::Map<Eigen::MatrixXcd>(ham_sq.data(), rows, cols);
//    double non_hermiticity = (ham_sq_map - ham_sq_map.adjoint()).cwiseAbs().sum() / static_cast<double>(ham_sq_map.size());
//    double sparcity        = static_cast<double>((ham_sq_map.array().cwiseAbs2() != 0.0).count()) / static_cast<double>(ham_sq_map.size());
//
//    if(non_hermiticity > 1e-12) throw std::runtime_error(fmt::format("multisite hamiltonian squared is not hermitian: {:.16f}", non_hermiticity));
//    if(non_hermiticity > 1e-14) tools::log->warn("multisite hamiltonian squared is slightly non-hermitian: {:.16f}", non_hermiticity);
//    if(ham_sq_map.hasNaN()) throw std::runtime_error("multisite hamiltonian squared has NaN's!");
//    tools::log->trace("multisite hamiltonian squared nonzeros: {:.8f} %", sparcity * 100);
//    tools::common::profile::t_hsq->toc();
//    return ham_sq;
//}

// Eigen::MatrixXcd tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_matrix(const class_model_finite &model,
//                                                                                                const class_edges_finite &edges) {
//    auto ham_tensor = tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_tensor(model, edges);
//    auto cols = ham_tensor.dimension(0) * ham_tensor.dimension(1) * ham_tensor.dimension(2);
//    auto rows = ham_tensor.dimension(3) * ham_tensor.dimension(4) * ham_tensor.dimension(5);
//    if(rows != cols) throw std::runtime_error(fmt::format("Hamiltonian tensor is not square: rows [{}] | cols [{}]",rows,cols));
//    return Eigen::Map<Eigen::MatrixXcd>(ham_tensor.data(), rows, cols).transpose().selfadjointView<Eigen::Lower>();
//}
//
// Eigen::MatrixXcd tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_matrix(const class_model_finite &model,
//                                                                                                        const class_edges_finite &edges) {
//    auto ham_squared_tensor = tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_squared_tensor(model, edges);
//    auto cols = ham_squared_tensor.dimension(0) * ham_squared_tensor.dimension(1) * ham_squared_tensor.dimension(2);
//    auto rows = ham_squared_tensor.dimension(3) * ham_squared_tensor.dimension(4) * ham_squared_tensor.dimension(5);
//    if(rows != cols) throw std::runtime_error(fmt::format("Hamiltonian squared tensor is not square: rows [{}] | cols [{}]",rows,cols));
//    return Eigen::Map<Eigen::MatrixXcd>(ham_squared_tensor.data(), rows, cols).transpose().selfadjointView<Eigen::Lower>();
//}

// Template definitions


//template<typename T>
//Eigen::Tensor<T, 6> tools::finite::opt::internal::get_multi_hamiltonian_tensor(const class_model_finite &model, const class_edges_finite &edges) {
//    static_assert(std::is_same<T, std::complex<double>>::value or std::is_same<T, double>::value, "Wrong type, expected double or complex double");
//    if constexpr(std::is_same<T, std::complex<double>>::value)
//        return tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_tensor(model, edges);
//    else if constexpr(std::is_same<T, double>::value)
//        return tools::finite::opt::internal::local_hamiltonians::get_multi_hamiltonian_tensor(model, edges).real();
//}

// Explicit instantiations
//template Eigen::Tensor<double, 6> tools::finite::opt::internal::get_multi_hamiltonian_tensor<double>(const class_model_finite &model,
//                                                                                                     const class_edges_finite &edges);
//template Eigen::Tensor<Scalar, 6> tools::finite::opt::internal::get_multi_hamiltonian_tensor<Scalar>(const class_model_finite &model,
//                                                                                                     const class_edges_finite &edges);

//template<typename T>
//Eigen::Tensor<T, 6> tools::finite::opt::internal::get_multi_hamiltonian_squared_tensor(const class_model_finite &model, const class_edges_finite &edges) {
//    static_assert(std::is_same<T, std::complex<double>>::value or std::is_same<T, double>::value, "Wrong type, expected double or complex double");
//    if constexpr(std::is_same<T, std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_squared_tensor(model, edges);
//    else if constexpr(std::is_same<T, double>::value)
//        return local_hamiltonians::get_multi_hamiltonian_squared_tensor(model, edges).real();
//}
// Explicit instantiations
//template Eigen::Tensor<double, 6> tools::finite::opt::internal::get_multi_hamiltonian_squared_tensor<double>(const class_model_finite &model,
//                                                                                                             const class_edges_finite &edges);
//template Eigen::Tensor<Scalar, 6> tools::finite::opt::internal::get_multi_hamiltonian_squared_tensor<Scalar>(const class_model_finite &model,
//                                                                                                             const class_edges_finite &edges);
//
//template<typename T>
//MatrixType<T> tools::finite::opt::internal::get_multi_hamiltonian_matrix(const class_model_finite &model, const class_edges_finite &edges) {
//    static_assert(std::is_same<T, std::complex<double>>::value or std::is_same<T, double>::value, "Wrong type, expected double or complex double");
//    if constexpr(std::is_same<T, std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_matrix(model, edges);
//    else if constexpr(std::is_same<T, double>::value)
//        return local_hamiltonians::get_multi_hamiltonian_matrix(model, edges).real();
//}
// Explicit instantiations

//
//template<typename T>
//MatrixType<T> tools::finite::opt::internal::get_multi_hamiltonian_squared_matrix(const class_model_finite &model, const class_edges_finite &edges) {
//    static_assert(std::is_same<T, std::complex<double>>::value or std::is_same<T, double>::value, "Wrong type, expected double or complex double");
//    if constexpr(std::is_same<T, std::complex<double>>::value) return local_hamiltonians::get_multi_hamiltonian_squared_matrix(model, edges);
//    else if constexpr(std::is_same<T, double>::value)
//        return local_hamiltonians::get_multi_hamiltonian_squared_matrix(model, edges).real();
//}
//// Explicit instantiations
//template MatrixType<double> tools::finite::opt::internal::get_multi_hamiltonian_squared_matrix<double>(const class_model_finite &model,
//                                                                                                       const class_edges_finite &edges);
//template MatrixType<Scalar> tools::finite::opt::internal::get_multi_hamiltonian_squared_matrix<Scalar>(const class_model_finite &model,
//                                                                                                       const class_edges_finite &edges);
//
//template<typename T>
//MatrixType<T> tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite &model, const class_edges_finite &edges,
//                                                                                              const std::vector<opt_tensor> &candidate_list) {
//    static_assert(std::is_same<T, std::complex<double>>::value or std::is_same<T, double>::value, "Wrong type, expected double or complex double");
//    if constexpr(std::is_same<T, std::complex<double>>::value)
//        return local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix_new(model, edges, candidate_list);
//    else if constexpr(std::is_same<T, double>::value)
//        return local_hamiltonians::get_multi_hamiltonian_squared_subspace_matrix_new(model, edges, candidate_list).real();
//}
//// Explicit instantiations
//template MatrixType<double>
//    tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new<double>(const class_model_finite &model, const class_edges_finite &edges,
//                                                                                            const std::vector<opt_tensor> &candidate_list);
//template MatrixType<Scalar>
//    tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new<Scalar>(const class_model_finite &model, const class_edges_finite &edges,
//                                                                                            const std::vector<opt_tensor> &candidate_list);
