#include <math/eig.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <config/enums.h>
template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace_full(const MatrixType<Scalar> & H_local,
                                                                                                         const TensorType<cplx, 3> &multisite_mps) {
    tools::log->trace("Finding subspace -- full");
    if(H_local.rows() != H_local.cols()) throw std::runtime_error(fmt::format("H_local is not square: rows [{}] | cols [{}]", H_local.rows(), H_local.cols()));
    auto t_opt_sub_eig = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_eig"]->tic_token();
    eig::solver solver;
    solver.eig<eig::Form::SYMM>(H_local.data(), H_local.rows(), eig::Vecs::ON, eig::Dephase::OFF);
    auto eigvals = eig::view::get_eigvals<double>(solver.result);
    auto eigvecs = eig::view::get_eigvecs<Scalar>(solver.result);
    tools::log->debug("Finished eigensolver -- reason: Full diagonalization");
    Eigen::Map<const Eigen::VectorXcd> multisite_vector(multisite_mps.data(), multisite_mps.size());
    Eigen::VectorXd                    overlaps = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
    int                                idx;
    double                             max_overlap    = overlaps.maxCoeff(&idx);
    double                             min_overlap    = overlaps.minCoeff();
    double                             sq_sum_overlap = overlaps.cwiseAbs2().sum();
    double                             subspace_error = 1.0 - sq_sum_overlap;
    long                               nev            = eigvecs.cols();
    reports::eigs_add_entry(nev, max_overlap, min_overlap, std::log10(subspace_error),
                            tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_ham"]->get_last_interval(),
                            tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_eig"]->get_last_interval(), 0, 1);
    return std::make_tuple(eigvecs, eigvals);
}

template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace_full(const MatrixType<cplx> &   H_local,
                                                                                                                  const TensorType<cplx, 3> &multisite_mps);

template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace_full(const MatrixType<real> &   H_local,
                                                                                                                  const TensorType<cplx, 3> &multisite_mps);
