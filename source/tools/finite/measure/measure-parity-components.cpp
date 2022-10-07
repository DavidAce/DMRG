#include "../measure.h"
#include "math/eig.h"
#include "math/linalg/matrix.h"
#include "math/linalg/tensor.h"
#include "math/tenx.h"
#include "qm/mpo.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
namespace tools::finite::measure {
    void parity_components(const StateFinite &state, const Eigen::Matrix2cd &paulimatrix) {
        auto t_spn                 = tid::tic_scope("spin", tid::level::highest);
        auto [mpo, L, R]           = qm::mpo::pauli_mpo(paulimatrix);
        auto                   mid = state.get_length<long>() / 2;
        Eigen::Tensor<cplx, 3> temp;
        tools::log->info("Labels: {}", state.get_labels());
        for(const auto &mps : state.mps_sites) {
            tools::common::contraction::contract_env_mps_mpo(temp, L, mps->get_M(), mpo);
            L = temp;
            if(mps->get_position<long>() == mid) {
                if(mps->get_label() == "AC") {
                    tools::log->info("Removing C");
                    temp = L.contract(tenx::asDiagonalInversed(mps->get_LC()), tenx::idx({0}, {0}))
                               .contract(tenx::asDiagonalInversed(mps->get_LC()), tenx::idx({0}, {0}));
                    L = temp;
                }
                if(mps->get_label() == "AC") {
                    tools::log->info("Removing L");
                    temp = L.contract(tenx::asDiagonalInversed(mps->get_L()), tenx::idx({0}, {0}))
                               .contract(tenx::asDiagonalInversed(mps->get_L()), tenx::idx({0}, {0}));
                    L = temp;
                }
                break;
            }
        }

        tools::log->info("Half spin product: \n{}\n", linalg::tensor::to_string(L.real(), 16, 20));
        eig::solver solver;
        solver.eig(L.data(), L.dimension(0));
        solver.config.tag = "zheevd";
        solver.eig<eig::Form::SYMM>(L.data(), L.dimension(0));
        tools::log->info("Half spin product eigenvalues: \n{}\n", linalg::matrix::to_string(eig::view::get_eigvals<eig::real>(solver.result), 16, 20));
    }
}
