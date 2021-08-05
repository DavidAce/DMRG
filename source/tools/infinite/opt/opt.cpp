#include "../opt.h"
#include <math/eig.h>
#include <math/eig/matvec/matvec_mpo.h>
#include <math/tenx.h>
#include <tensors/edges/EdgesInfinite.h>
#include <tensors/model/ModelInfinite.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateInfinite.h>
#include <tensors/TensorsInfinite.h>
#include <tid/tid.h>

namespace tools::infinite::opt {
    Eigen::Tensor<cplx, 3> find_ground_state(const TensorsInfinite &state, OptRitz ritz) {
        return tools::infinite::opt::find_ground_state(state, enum2sv(ritz));
    }

    Eigen::Tensor<cplx, 3> find_ground_state(const TensorsInfinite &tensors, std::string_view ritzstring) {
        tools::log->trace("Starting ground state optimization");
        auto t_eig = tid::tic_scope("eig");

        eig::Ritz ritz = eig::stringToRitz(ritzstring);

        auto        shape_mps = tensors.state->dimensions();
        auto        shape_mpo = tensors.model->dimensions();
        const auto &mpo       = tensors.model->get_2site_mpo_AB();
        const auto &env       = tensors.edges->get_ene_blk();

        MatVecMPO<cplx> matrix(env.L.data(), env.R.data(), mpo.data(), shape_mps, shape_mpo);
        eig::solver     solver;
        solver.config.maxNev  = 1;
        solver.config.maxNcv  = static_cast<eig::size_type>(settings::precision::eig_default_ncv);
        solver.config.tol     = settings::precision::eig_tolerance;
        solver.config.maxIter = 10000;
        solver.eigs(matrix, -1, -1, ritz, eig::Form::SYMM, eig::Side::R, 1.0, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF);
        return eig::view::get_eigvec<cplx>(solver.result, shape_mps);
    }

    //============================================================================//
    // Do unitary evolution on an MPS
    //============================================================================//

    Eigen::Tensor<cplx, 3> time_evolve_state(const StateInfinite &state, const Eigen::Tensor<cplx, 2> &U)
    /*!
    @verbatim
      1--[ mps ]--2
            |
            0
                             1--[ mps ]--2
            0         --->         |
            |                      0
          [ U ]
            |
            1
    @endverbatim
    */
    {
        return U.contract(state.get_2site_mps(), tenx::idx({0}, {0}));
    }

}
