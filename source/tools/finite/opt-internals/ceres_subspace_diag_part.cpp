#include <config/nmspc_settings.h>
#include <math/arpack_extra/matrix_product_stl.h>
#include <math/class_eigsolver.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt.h>
#include <vector>

std::vector<int> tools::finite::opt::internal::generate_size_list(int shape) {
    int max_nev;
    if(shape <= 512) max_nev = shape / 2;
    else if(shape > 512 and shape <= 1024)
        max_nev = shape / 4;
    // should do full diag
    else if(shape > 1024 and shape <= 2048)
        max_nev = shape / 4;
    // should do full diag
    else if(shape > 2048 and shape <= 4096)
        max_nev = 64;
    else if(shape > 4096 and shape <= 8192)
        max_nev = 32;
    else
        max_nev = 16;

    int min_nev = std::min(std::min(8, (int) shape), max_nev);

    std::vector<int> nev_list = {min_nev};
    int              tmp_nev  = min_nev;
    while(tmp_nev < max_nev) {
        tmp_nev = std::min(4 * tmp_nev, max_nev);
        nev_list.push_back(tmp_nev);
    }
    return nev_list;
}

template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
    tools::finite::opt::internal::subspace::find_subspace_part(const MatrixType<Scalar> &H_local, const TensorType<cplx, 3> &multisite_tensor,
                                                               double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace) {
    using namespace eigutils::eigSetting;
    using namespace tools::common::profile;
    tools::log->trace("Finding subspace -- partial");
    // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
    StlMatrixProduct<Scalar> hamiltonian(H_local.data(), static_cast<int>(H_local.rows()), true, Form::SYMMETRIC, Side::R);
    hamiltonian.set_shift(energy_target);
    hamiltonian.FactorOP();
    double time_lu  = hamiltonian.t_factorOp.get_last_time_interval();
    double time_ham = tools::common::profile::t_ham->get_last_time_interval();
    tools::common::profile::t_eig->toc();

    class_eigsolver solver;
    solver.solverConf.eigThreshold            = settings::precision::eig_threshold;
    std::string                        reason = "exhausted";
    Eigen::VectorXd                    eigvals;
    Eigen::MatrixXcd                   eigvecs;
    Eigen::Map<const Eigen::VectorXcd> multisite_vector(multisite_tensor.data(), multisite_tensor.size());
    for(auto nev : generate_size_list(static_cast<int>(multisite_tensor.size()))) {
        t_eig->tic();
        solver.eigs_stl(hamiltonian, nev, -1, energy_target, Form::SYMMETRIC, Ritz::LM, Side::R, true, false);
        t_eig->toc();

        eigvals = Eigen::Map<const Eigen::VectorXd>(solver.solution.get_eigvals<Form::SYMMETRIC>().data(), solver.solution.meta.cols);
        if constexpr(std::is_same<std::complex<double>, Scalar>::value) {
            eigvecs = Eigen::Map<const Eigen::MatrixXcd>(solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(), solver.solution.meta.rows,
                                                         solver.solution.meta.cols);
        } else {
            eigvecs = Eigen::Map<const Eigen::MatrixXd>(solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(), solver.solution.meta.rows,
                                                        solver.solution.meta.cols);
        }

        Eigen::VectorXd overlaps       = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
        double          max_overlap    = overlaps.maxCoeff();
        double          min_overlap    = overlaps.minCoeff();
        double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
        double          subspace_error = 1.0 - sq_sum_overlap;
        reports::eigs_log.emplace_back(nev, max_overlap, min_overlap, std::log10(subspace_error), time_ham, t_eig->get_last_time_interval(), time_lu);
        time_lu  = 0;
        time_ham = 0;
        if(max_overlap > 1.0 + 1e-6) throw std::runtime_error("max_overlap larger than one : " + std::to_string(max_overlap));
        if(sq_sum_overlap > 1.0 + 1e-6) throw std::runtime_error("eps larger than one : " + std::to_string(sq_sum_overlap));
        if(min_overlap < 0.0) throw std::runtime_error("min_overlap smaller than zero: " + std::to_string(min_overlap));
        if(subspace_error < subspace_error_threshold) {
            reason = "subspace error is low enough";
            break;
        }
        if(optSpace == OptSpace::SUBSPACE_AND_DIRECT and subspace_error < 1e-3) {
            reason = "subspace error sufficient for SUBSPACE_AND_DIRECT";
            break;
        }
        if(optSpace == OptSpace::SUBSPACE_ONLY and optMode == OptMode::OVERLAP and max_overlap >= 1.0 / std::sqrt(2.0)) {
            reason = "Overlap sufficient for OVERLAP and SUBSPACE_ONLY";
            break;
        }
    }
    tools::log->debug("Finished partial eigensolver -- reason: {}", reason);
    tools::common::profile::t_eig->toc();
    return std::make_tuple(eigvecs, eigvals);
}


template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
tools::finite::opt::internal::subspace::find_subspace_part(const MatrixType<real> &H_local, const TensorType<cplx,3> &multisite_tensor, double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
tools::finite::opt::internal::subspace::find_subspace_part(const MatrixType<cplx> &H_local, const TensorType<cplx,3> &multisite_tensor, double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

