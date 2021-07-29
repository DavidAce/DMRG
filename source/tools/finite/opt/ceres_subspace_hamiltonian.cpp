#include <math/tenx.h>

// -- (textra first)

#include <tensors/edges/EdgesFinite.h>
#include <tensors/model/ModelFinite.h>
#include <tid/tid.h>
#include <tools/common/contraction.h>
#include <tools/common/log.h>
#include <tools/finite/opt/opt-internal.h>
#include <tools/finite/opt_mps.h>

using Scalar = std::complex<double>;
using real   = double;
using cplx   = std::complex<double>;
template<typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar>
void shift_mpo_energy(Eigen::Tensor<Scalar, 4> &mpo, double energy_shift) {
    tools::log->debug("Setting shift: {:.16f}", energy_shift);

    // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
    // rank2 matrix, where each element is also a matrix with the size
    // determined by the last 2 indices kl.
    // When we shift an MPO, all we do is subtract a diagonal matrix from
    // the botton left corner of the ij-matrix.
    // Setup extents and handy objects
    auto                                       shape = mpo.dimensions();
    std::array<long, 4>                        offset4{shape[0] - 1, 0, 0, 0};
    std::array<long, 4>                        extent4{1, 1, shape[2], shape[3]};
    std::array<long, 2>                        extent2{shape[2], shape[3]};
    MatrixType<Scalar>                         sigma_Id = energy_shift * MatrixType<Scalar>::Identity(extent2[0], extent2[1]);
    Eigen::TensorMap<Eigen::Tensor<Scalar, 2>> sigma_Id_map(sigma_Id.data(), sigma_Id.rows(), sigma_Id.cols());
    mpo.slice(offset4, extent4).reshape(extent2) -= sigma_Id_map;
}

template<typename T>
MatrixType<T> tools::finite::opt::internal::get_multisite_hamiltonian_matrix(const ModelFinite &model, const EdgesFinite &edges, double energy_shift) {
    auto        t_ham = tid::tic_scope("ham");
    const auto &env   = edges.get_multisite_ene_blk();
    auto        mpo   = model.get_multisite_mpo();
    if(energy_shift != 0.0) shift_mpo_energy(mpo, energy_shift);
    tools::log->trace("Contracting multisite hamiltonian");

    long dim0 = mpo.dimension(2);
    long dim1 = env.L.dimension(0);
    long dim2 = env.R.dimension(0);

    Eigen::Tensor<std::complex<double>, 6> ham(dim0, dim1, dim2, dim0, dim1, dim2);
    ham.device(tenx::omp::getDevice()) = env.L.contract(mpo, tenx::idx({2}, {0})).contract(env.R, tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});

    auto cols = ham.dimension(0) * ham.dimension(1) * ham.dimension(2);
    auto rows = ham.dimension(3) * ham.dimension(4) * ham.dimension(5);
    long size = dim0 * dim1 * dim2;
    if(rows != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian rows (dim0*dim1*dim2) and size: {} != {}", cols, size));
    if(cols != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian cols (dim3*dim4*dim5) and size: {} != {}", rows, size));

    auto   ham_map         = Eigen::Map<Eigen::MatrixXcd>(ham.data(), rows, cols);
    double non_hermiticity = (ham_map - ham_map.adjoint()).cwiseAbs().sum() / static_cast<double>(ham_map.size());
    double sparcity        = static_cast<double>((ham_map.array().cwiseAbs2() != 0.0).count()) / static_cast<double>(ham_map.size());

    if(non_hermiticity > 1e-8) throw std::runtime_error(fmt::format("multisite hamiltonian is not hermitian: {:.16f}", non_hermiticity));
    if(non_hermiticity > 1e-12) tools::log->warn("multisite hamiltonian is slightly non-hermitian: {:.16f}", non_hermiticity);
    if(ham_map.hasNaN()) throw std::runtime_error("multisite hamiltonian has NaN's!");
    tools::log->trace("multisite hamiltonian nonzeros: {:.8f} %", sparcity * 100);
    if constexpr(std::is_same_v<T, double>) {
        return Eigen::Map<Eigen::MatrixXcd>(ham.data(), rows, cols).real().transpose().selfadjointView<Eigen::Lower>();
    } else {
        return Eigen::Map<Eigen::MatrixXcd>(ham.data(), rows, cols).transpose().selfadjointView<Eigen::Lower>();
    }
}
template MatrixType<real> tools::finite::opt::internal::get_multisite_hamiltonian_matrix(const ModelFinite &model, const EdgesFinite &edges,
                                                                                         double energy_shift);
template MatrixType<cplx> tools::finite::opt::internal::get_multisite_hamiltonian_matrix(const ModelFinite &model, const EdgesFinite &edges,
                                                                                         double energy_shift);

template<typename T>
MatrixType<T> tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix(const ModelFinite &model, const EdgesFinite &edges) {
    auto        t_opt_sub_ham = tid::tic_scope("ham²");
    const auto &mpo2          = model.get_multisite_mpo_squared();
    const auto &env2          = edges.get_multisite_var_blk();
    tools::log->trace("Contracting multisite hamiltonian");

    long dim0 = mpo2.dimension(2);
    long dim1 = env2.L.dimension(0);
    long dim2 = env2.R.dimension(0);

    Eigen::Tensor<std::complex<double>, 6> ham2(dim0, dim1, dim2, dim0, dim1, dim2);
    ham2.device(tenx::omp::getDevice()) =
        env2.L.contract(mpo2, tenx::idx({2}, {0})).contract(env2.R, tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});

    auto cols = ham2.dimension(0) * ham2.dimension(1) * ham2.dimension(2);
    auto rows = ham2.dimension(3) * ham2.dimension(4) * ham2.dimension(5);
    long size = dim0 * dim1 * dim2;
    if(rows != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian rows (dim0*dim1*dim2) and size: {} != {}", cols, size));
    if(cols != size) throw std::runtime_error(fmt::format("Mismatch in multisite hamiltonian cols (dim3*dim4*dim5) and size: {} != {}", rows, size));

    auto   ham2_map        = Eigen::Map<Eigen::MatrixXcd>(ham2.data(), rows, cols);
    double non_hermiticity = (ham2_map - ham2_map.adjoint()).cwiseAbs().sum() / static_cast<double>(ham2_map.size());
    double sparcity        = static_cast<double>((ham2_map.array().cwiseAbs2() != 0.0).count()) / static_cast<double>(ham2_map.size());

    if(non_hermiticity > 1e-8) throw std::runtime_error(fmt::format("multisite hamiltonian squared is not hermitian: {:.16f}", non_hermiticity));
    if(non_hermiticity > 1e-12) tools::log->warn("multisite hamiltonian squared is slightly non-hermitian: {:.16f}", non_hermiticity);
    if(ham2_map.hasNaN()) throw std::runtime_error("multisite hamiltonian squared has NaN's!");
    tools::log->trace("multisite hamiltonian squared nonzeros: {:.8f} %", sparcity * 100);
    if constexpr(std::is_same_v<T, double>) {
        return Eigen::Map<Eigen::MatrixXcd>(ham2.data(), rows, cols).real().transpose().selfadjointView<Eigen::Lower>();
    } else {
        return Eigen::Map<Eigen::MatrixXcd>(ham2.data(), rows, cols).transpose().selfadjointView<Eigen::Lower>();
    }
}
template MatrixType<real> tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix(const ModelFinite &model, const EdgesFinite &edges);
template MatrixType<cplx> tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix(const ModelFinite &model, const EdgesFinite &edges);

template<typename T>
MatrixType<T> tools::finite::opt::internal::get_multisite_hamiltonian_squared_subspace_matrix(const ModelFinite &model, const EdgesFinite &edges,
                                                                                              const std::vector<opt_mps> &candidate_list, double energy_shift) {
    // First, make sure every candidate is actually a basis vector, otherwise this computation would turn difficult if we have to skip rows and columns
    auto t_ham = tid::tic_scope("ham²_sub");
    for(const auto &candidate : candidate_list)
        if(not candidate.is_basis_vector)
            throw std::runtime_error("One candidate is not a basis vector. When constructing a hamiltonian subspace matrix, make sure the candidates are all "
                                     "eigenvectors/basis vectors");

    const auto &env2 = edges.get_multisite_var_blk();
    auto        mpo2 = model.get_multisite_mpo_squared();
    if(energy_shift != 0.0 and not model.is_compressed_mpo_squared()) shift_mpo_energy(mpo2, energy_shift);

    tools::log->trace("Contracting subspace hamiltonian squared new");
    long dim0   = mpo2.dimension(2);
    long dim1   = env2.L.dimension(0);
    long dim2   = env2.R.dimension(0);
    long eignum = static_cast<long>(candidate_list.size()); // Number of eigenvectors

    Eigen::Tensor<Scalar, 0> H2_ij;
    Eigen::Tensor<Scalar, 3> H2_mps(dim0, dim1, dim2); // The local hamiltonian multiplied by mps at column j.
    Eigen::MatrixXcd         H2_sub(eignum, eignum);   // The local hamiltonian projected to the subspace (spanned by candidate_list)

    for(auto col = 0; col < eignum; col++) {
        const auto &mps_j = std::next(candidate_list.begin(), col)->get_tensor();
        tools::common::contraction::matrix_vector_product(H2_mps, mps_j, mpo2, env2.L, env2.R);
        for(auto row = col; row < eignum; row++) {
            const auto &mps_i                    = std::next(candidate_list.begin(), row)->get_tensor();
            H2_ij.device(tenx::omp::getDevice()) = mps_i.conjugate().contract(H2_mps, tenx::idx({0, 1, 2}, {0, 1, 2}));
            H2_sub(row, col)                     = H2_ij(0);
            H2_sub(col, row)                     = std::conj(H2_ij(0));
        }
    }
    if constexpr(std::is_same_v<T, double>)
        return H2_sub.real();
    else
        return H2_sub;
}

// Explicit instantiations
template MatrixType<real> tools::finite::opt::internal::get_multisite_hamiltonian_squared_subspace_matrix(const ModelFinite &model, const EdgesFinite &edges,
                                                                                                          const std::vector<opt_mps> &candidate_list,
                                                                                                          double                      energy_shift);
template MatrixType<cplx> tools::finite::opt::internal::get_multisite_hamiltonian_squared_subspace_matrix(const ModelFinite &model, const EdgesFinite &edges,
                                                                                                          const std::vector<opt_mps> &candidate_list,
                                                                                                          double                      energy_shift);
