#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include <math/class_eigsolver.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt.h>



template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace_full(const MatrixType<Scalar> &H_local,
                                                                                                         const TensorType<cplx, 3> &   multisite_tensor) {
    using namespace eigutils::eigSetting;
    using namespace tools::common::profile;
    tools::log->trace("Finding subspace -- full");
    if(H_local.rows() != H_local.cols()) throw std::runtime_error(fmt::format("H_local is not square: rows [{}] | cols [{}]", H_local.rows(), H_local.cols()));
    Eigen::VectorXd  eigvals;
    Eigen::MatrixXcd eigvecs;
    class_eigsolver  solver;

    t_eig->tic();
    if constexpr(!std::is_same<Scalar, double>::value) {
        solver.eig<Type::CPLX, Form::SYMMETRIC>(H_local, true, false);
        eigvals = Eigen::Map<const Eigen::VectorXd>(solver.solution.get_eigvals<Form::SYMMETRIC>().data(), solver.solution.meta.cols);
        eigvecs = Eigen::Map<const Eigen::MatrixXcd>(solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(), solver.solution.meta.rows,
                                                     solver.solution.meta.cols);
    } else {
        solver.eig<Type::REAL, Form::SYMMETRIC>(H_local, true, false);
        eigvals = Eigen::Map<const Eigen::VectorXd>(solver.solution.get_eigvals<Form::SYMMETRIC>().data(), solver.solution.meta.cols);
        eigvecs = Eigen::Map<const Eigen::MatrixXd>(solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(), solver.solution.meta.rows,
                                                    solver.solution.meta.cols);
    }
    t_eig->toc();
    tools::log->debug("Finished eigensolver -- reason: Full diagonalization");

    Eigen::Map<const Eigen::VectorXcd> theta_vec(multisite_tensor.data(), multisite_tensor.size());
    Eigen::VectorXd                    overlaps = (theta_vec.adjoint() * eigvecs).cwiseAbs().real();
    int                                idx;
    double                             max_overlap    = overlaps.maxCoeff(&idx);
    double                             min_overlap    = overlaps.minCoeff();
    double                             sq_sum_overlap = overlaps.cwiseAbs2().sum();
    double                             subspace_error = 1.0 - sq_sum_overlap;
    int                                nev            = static_cast<int>(eigvecs.cols());
    reports::eigs_log.emplace_back(nev, max_overlap, min_overlap, std::log10(subspace_error), t_ham->get_last_time_interval(), t_eig->get_last_time_interval(),
                                   0);
    return std::make_tuple(eigvecs, eigvals);
}



template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace_full(const MatrixType<cplx> &H_local,
                                                                                                                  const TensorType<cplx, 3> & multisite_tensor);

template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace_full(const MatrixType<real> &H_local,
                                                                                                                  const TensorType<cplx, 3> & multisite_tensor);

