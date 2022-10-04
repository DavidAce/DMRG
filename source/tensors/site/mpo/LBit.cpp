#include "LBit.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/linalg/tensor.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/spin.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>

namespace settings {
    inline constexpr bool debug_lbit = false;
}

std::string format_mpo(const Eigen::Tensor<std::complex<double>, 4> &mpo, double J1_rand, const std::vector<double> &J2_rand, double J3_rand) {
    std::string str;
    using namespace qm::spin::half;
    auto ext = std::array<long, 4>{1, 1, 2, 2};
    auto sh2 = std::array<long, 2>{2, 2};
    for(long i = 0; i < mpo.dimension(0); i++) {
        for(long j = 0; j < mpo.dimension(1); j++) {
            auto                                   off = std::array<long, 4>{i, j, 0, 0};
            Eigen::Tensor<std::complex<double>, 2> mat = mpo.slice(off, ext).reshape(sh2);
            if(tenx::MatrixMap(mat) == sz)
                str.append("z   ");
            else if(tenx::MatrixMap(mat) == id)
                str.append("i   ");
            else if(tenx::isZero(mat))
                str.append("0   ");
            else {
                if(tenx::MatrixMap(mat).isApprox(J1_rand * sz)) str.append(fmt::format("J1z "));
                if(tenx::MatrixMap(mat).isApprox(J3_rand * sz)) str.append(fmt::format("J3z "));
                for(const auto &[r, J2] : iter::enumerate(J2_rand)) {
                    if(tenx::MatrixMap(mat).isApprox(J2 * sz)) str.append(fmt::format("J2{} ", r));
                }
            }
        }
        str.append("\n");
    }
    return str;
}

LBit::LBit(ModelType model_type_, size_t position_) : MpoSite(model_type_, position_) {
    h5tb.param.J1_mean  = settings::model::lbit::J1_mean;
    h5tb.param.J2_mean  = settings::model::lbit::J2_mean;
    h5tb.param.J3_mean  = settings::model::lbit::J3_mean;
    h5tb.param.J1_wdth  = settings::model::lbit::J1_wdth;
    h5tb.param.J2_wdth  = settings::model::lbit::J2_wdth;
    h5tb.param.J3_wdth  = settings::model::lbit::J3_wdth;
    h5tb.param.J2_xcls  = settings::model::lbit::J2_xcls;
    h5tb.param.J2_span  = settings::model::lbit::J2_span; // Can be 0. If -1ul = MAX then this means we want the cutoff to be full system range
    h5tb.param.f_mixer  = settings::model::lbit::f_mixer;
    h5tb.param.u_layer  = settings::model::lbit::u_layer;
    h5tb.param.spin_dim = settings::model::lbit::spin_dim;

    // Adjust J2_span, it doesn't make sense to have it larger than the system size anyway, so we use a cutoff
    h5tb.param.J2_ctof      = std::min(h5tb.param.J2_span, settings::model::model_size - 1); // Range unsigned long
    h5tb.param.distribution = settings::model::lbit::distribution;
    extent4                 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2                 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    all_mpo_parameters_have_been_set =
        false; // There are no full lattice parameters on lbit, but we set it to false so we remember to call randomize on all sites in every model type
}

// double LBit::get_field() const { return h5tb.param.J1_rand; }
// double LBit::get_coupling() const { return h5tb.param.J2_rand + h5tb.param.J3_rand; }
void LBit::print_parameter_names() const { h5tb.print_parameter_names(); }
void LBit::print_parameter_values() const { h5tb.print_parameter_values(); }

void LBit::set_parameters(TableMap &parameters) {
    h5tb.param.J1_rand      = std::any_cast<double>(parameters["J1_rand"]);
    h5tb.param.J2_rand      = std::any_cast<h5pp::varr_t<double>>(parameters["J2_rand"]);
    h5tb.param.J3_rand      = std::any_cast<double>(parameters["J3_rand"]);
    h5tb.param.J1_wdth      = std::any_cast<double>(parameters["J1_wdth"]);
    h5tb.param.J2_wdth      = std::any_cast<double>(parameters["J2_wdth"]);
    h5tb.param.J3_wdth      = std::any_cast<double>(parameters["J3_wdth"]);
    h5tb.param.J2_xcls      = std::any_cast<double>(parameters["J2_xcls"]);
    h5tb.param.J2_span      = std::any_cast<size_t>(parameters["J2_span"]);
    h5tb.param.f_mixer      = std::any_cast<double>(parameters["f_mixer"]);
    h5tb.param.u_layer      = std::any_cast<size_t>(parameters["u_layer"]);
    h5tb.param.spin_dim     = std::any_cast<long>(parameters["spin_dim"]);
    h5tb.param.distribution = std::any_cast<h5pp::vstr_t>(parameters["distribution"]);
    // Adjust J2_span, it doesn't make sense to have it larger than the system size anyway, so we use a cutoff
    h5tb.param.J2_ctof               = std::min(h5tb.param.J2_span, settings::model::model_size - 1); // Range unsigned long
    all_mpo_parameters_have_been_set = true;
}

LBit::TableMap LBit::get_parameters() const {
    /* clang-format off */
    TableMap parameters;
    parameters["J1_rand"]       = h5tb.param.J1_rand;
    parameters["J2_rand"]       = h5tb.param.J2_rand;
    parameters["J3_rand"]       = h5tb.param.J3_rand;
    parameters["J1_wdth"]       = h5tb.param.J1_wdth;
    parameters["J2_wdth"]       = h5tb.param.J2_wdth;
    parameters["J3_wdth"]       = h5tb.param.J3_wdth;
    parameters["J2_xcls"]       = h5tb.param.J2_xcls;
    parameters["J2_span"]       = h5tb.param.J2_span;
    parameters["J2_ctof"]       = h5tb.param.J2_ctof;
    parameters["f_mixer"]       = h5tb.param.f_mixer;
    parameters["u_layer"]       = h5tb.param.u_layer;
    parameters["spin_dim"]      = h5tb.param.spin_dim;
    parameters["distribution"]  = h5tb.param.distribution;
    return parameters;
    /* clang-format on */
}

void LBit::build_mpo()

/*
      Builds the MPO for the l-bit Hamiltonian as a rank-4 tensor.

      H =   Σ_i  J1(i)   * n_{i}
          + Σ_ij J2(i,j) * n_{i} * n_{i+j}
          + Σ_i  J3(i)   * n_{i} * n_{i+1} * n_{i+2}


      MPO example with F = 10,  R = 8:
                             0     1     2     3     4     5     6     7     8    9    F
           2            |    I     .     .     .     .     .     .     .     .    .    . | 0
           |            |    n     .     .     .     .     .     .     .     .    .    . | 1
       0---M---1    =   |    .     I     .     .     .     .     .     .     .    .    . | 2
           |            |    .     .     I     .     .     .     .     .     .    .    . | 3
           3            |    .     .     .     I     .     .     .     .     .    .    . | 4
                        |    .     .     .     .     I     .     .     .     .    .    . | 5
                        |    .     .     .     .     .     I     .     .     .    .    . | 6
                        |    .     .     .     .     .     .     I     .     .    .    . | 7
                        |    .     .     .     .     .     .     .     I     .    .    . | 8
                        |    .     n     .     .     .     .     .     .     .    .    . | 9
                        | J1*n J21*n J22*n J23*n J24*n J25*n J26*n J27*n J28*n J3*n    I | F

       where
        *   F is the linear size of the MPO
        *   R is the the max range of interactions, i.e site i can interact with i+8
        *	i,j are site indices
        *   n = 0.5 * (1 + σ^z),
        *   σ^z is the diagonal 2x2 pauli matrix
        *   I is the 2x2 identity matrix
        *   J1,J2? and J3 are random 1,2 and 3-body couplings
        *   Jij = J21, J22... couples sites i,j at |i-j| <= 5
        *   Jij = J2_rand(i,j) = exp(-|i-j|/J2_xcls) * Random(-w+m,w+m) , where
                * J2_decay is the exponential decay rate (with respect to distance) of the 2-body interaction
                * w is the box width of the uniform distribution
                * m=mean is a constant offset of the distribution

       The MPO is built from the following finite-state-machine ([k] are matrix indices and the machine gives us an adjacency matrix)

       I==[0]-------------J1*n(i)--------------[F]==I
           |                                    |
           |---n(i)---[1]----------J21*n(i+1)---| if R >= 1
                       |                        |
                       |-n(i+1)-[F-1]-J3*n(i+2)-|
                       |                        |
                       |--I--[2]---J22*n(i+2)---| if R >= 2
                              I                 |
                             [3]---J23*n(i+3)---| if R >= 3
                              I                 |
                             [4]---J24*n(i+4)---| ...
                              I                 |
                             [5]---J25*n(i+5)---|
                              I                 |
                             [6]---J26*n(i+6)---|
                              I                 |
                             [7]---J27*n(i+7)---|
                              I                 |
                             [8]---J28*n(i+8)---|
                             .
                             .
                             .

       The final step of the state machine, "F" depends on the longest pairwise range "R = max|i-j|" that we have in the system
            F = 1(for 1-body) + 1(for 3-body) + R (for all the 2-body)


        We see that
        - dims(M) = {F+1,F+1}
        - dims(J2) = {R+1}
        - nonzero elements are:
            - M[0,0] = M[F,F] = I
            - M[1,0] = n
            - M[i,i-1] = I for i in range(2,R+1) <--- "range" excludes the last element
            - M[F-1,1] = n
            - M[F,0] =  J1*n
            - M[F,i] = J2[i] for i in range(1,R+1) <--- "range" excludes the last element
            - M[F,F-1] = J3*n

  */

{
    using namespace qm::spin::half;
    tools::log->debug("mpo({}): building lbit mpo", get_position());
    if(not all_mpo_parameters_have_been_set)
        throw except::runtime_error("mpo({}): can't build mpo: full lattice parameters haven't been set yet.", get_position());
    if(h5tb.param.J2_rand.empty()) throw except::logic_error("expected non-empty J2_rand");
    if(h5tb.param.J2_rand.length() > h5tb.param.J2_ctof + 1)
        throw except::logic_error("expected J2_rand.length()({}) <= 1+J2_ctof ({})", h5tb.param.J2_rand.length(), h5tb.param.J2_ctof);

//    Eigen::Tensor<cplx, 2> n = tenx::TensorCast(0.5 * (id + sz));
#pragma message "using sz instead of number operator"
    Eigen::Tensor<cplx, 2> n = tenx::TensorMap(sz);
    Eigen::Tensor<cplx, 2> I = tenx::TensorMap(id);
    long                   R = static_cast<long>(h5tb.param.J2_ctof);
    long                   F = R + 2l;
    mpo_internal.resize(F + 1, F + 1, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();

    mpo_internal.slice(tenx::array4{0, 0, 0, 0}, extent4).reshape(extent2) = I;
    mpo_internal.slice(tenx::array4{F, F, 0, 0}, extent4).reshape(extent2) = I;
    mpo_internal.slice(tenx::array4{1, 0, 0, 0}, extent4).reshape(extent2) = n;

    if(R >= 2)
        for(const auto &i : num::range<long>(2, R + 1)) { mpo_internal.slice(tenx::array4{i, i - 1, 0, 0}, extent4).reshape(extent2) = I; }

    mpo_internal.slice(tenx::array4{F - 1, 1, 0, 0}, extent4).reshape(extent2) = n;
    mpo_internal.slice(tenx::array4{F, 0, 0, 0}, extent4).reshape(extent2)     = h5tb.param.J1_rand * n - e_shift * I;

    if(R >= 1)
        for(const auto &i : num::range<long>(1, h5tb.param.J2_rand.length())) {
            mpo_internal.slice(tenx::array4{F, i, 0, 0}, extent4).reshape(extent2) = h5tb.param.J2_rand[static_cast<size_t>(i)] * n;
        }
    mpo_internal.slice(tenx::array4{F, F - 1, 0, 0}, extent4).reshape(extent2) = h5tb.param.J3_rand * n;

    if(tenx::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw except::runtime_error("mpo({}): found nan", get_position());
    }
    unique_id = std::nullopt;
}

void LBit::randomize_hamiltonian() {
    // J2(i,j) = J2_exp(r) * Random_ij(J2_mean, J2_wdth), for i < j,
    // where
    //    * J2_exp(r) = exp(-(r-1)/J2_xcls)
    //    * J2_xcls: characteristic length scale for decay of pairwise interactions
    //    * r = r(i,j) = |i-j|
    //    * Random_ij(J2_mean, J2_wdth) are drawn randomly for each i,j, according to some distribution.
    //    * The "-1" in the exponent disables the exponential suppression on nearest neighbors,
    //      i.e.  exponential decay starts after 1 site, and so J2_wdth sets the size of nearest neighbor interaction.
    //
    // For normal and lognormal we can instead put the decay in the standard deviation
    //      J2(i,j) = N(J2_mean, J2_exp(r))
    // This definition is the same as Eq. 7: https://link.aps.org/doi/10.1103/PhysRevB.97.214202
    // Note that the saturation time is predictably:
    //      tmax ~ [J2_wdth * exp(-(L/2 - 1)/J2_xcls)]^-1
    // where 2/pi comes from using the half-normal distribution |N(...)| (see wiki).
    // We take the absolute value here to get the order of magnitude of the interactions.

    using namespace settings::model::lbit;
    auto J2_rnd = std::vector<double>();
    auto J2_exp = std::vector<double>();
    for(size_t r = 0; r < h5tb.param.J2_ctof + 1; ++r) {
        if(r == 0) {
            // J2 does not describe self-interaction
            J2_exp.emplace_back(0);
            continue;
        }
        if(r + get_position() >= settings::model::model_size) break; // No more sites to interact with
        J2_exp.emplace_back(std::exp(-(static_cast<double>(r) - 1) / settings::model::lbit::J2_xcls));
    }
    if(h5tb.param.distribution == "normal") {
        h5tb.param.J1_rand = rnd::normal(J1_mean, J1_wdth);
        h5tb.param.J3_rand = rnd::normal(J3_mean, J3_wdth);
        for(const auto &w : J2_exp) J2_rnd.emplace_back(rnd::normal(J2_mean, w));
    } else if(h5tb.param.distribution == "lognormal") {
        h5tb.param.J1_rand = rnd::log_normal(J1_mean, J1_wdth);
        h5tb.param.J3_rand = rnd::log_normal(J3_mean, J3_wdth);
        for(const auto &w : J2_exp) J2_rnd.emplace_back(rnd::log_normal(J2_mean, w));
    } else if(h5tb.param.distribution == "uniform") {
        h5tb.param.J1_rand = rnd::uniform_double_box(J1_mean - J1_wdth / 2.0, J1_mean + J1_wdth / 2.0);
        h5tb.param.J3_rand = rnd::uniform_double_box(J3_mean - J3_wdth / 2.0, J3_mean + J3_wdth / 2.0);
        for(const auto &w : J2_exp) J2_rnd.emplace_back(w * rnd::uniform_double_box(J2_mean - J2_wdth / 2.0, J2_mean + J2_wdth / 2.0));
    } else if(h5tb.param.distribution == "constant") {
        h5tb.param.J1_rand = settings::model::lbit::J1_mean;
        h5tb.param.J3_rand = settings::model::lbit::J3_mean;
        for(const auto &w : J2_exp) J2_rnd.emplace_back(w * settings::model::lbit::J2_mean);
    } else {
        throw except::runtime_error("wrong distribution [{}]: expected one of normal | lognormal | uniform | constant", h5tb.param.distribution);
    }
    h5tb.param.J2_rand               = J2_rnd;
    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
}

Eigen::Tensor<MpoSite::cplx, 4> LBit::MPO_nbody_view(std::optional<std::vector<size_t>> nbody, std::optional<std::vector<size_t>> skip) const {
    // This function returns a view of the MPO including only n-body terms.
    // For instance, if nbody_terms == {2,3}, this would include 2-body and 3-body but exclude 1-body on-site terms.

    // Meanwhile, skip is a list of mpo positions that we should ignore. This lets us disable interactions with certain sites.
    // For instance if skip == {2}, then interaction terms such as J[0,2], J[1,2] and J[2,3] are set to zero.

    if(not nbody) return MPO();
    auto   Rul                        = h5tb.param.J2_ctof;                    // Range unsigned long (with "full range": R = L-1)
    long   R                          = static_cast<long>(h5tb.param.J2_ctof); // Range
    long   F                          = R + 2l;                                // Last index of mpo
    size_t pos                        = get_position();
    auto   J2_range                   = num::range<size_t>(1, R + 1); // +1 so that R itself is included
    double J1_on                      = 0.0;
    double J2_on                      = 0.0;
    double J3_on                      = 0.0;
    bool   adjust_for_double_counting = false;
    // Toggle on requested nbody terms
    for(const auto &n : nbody.value()) {
        if(n == 0) adjust_for_double_counting = true;
        if(n == 1) J1_on = 1.0;
        if(n == 2) J2_on = 1.0;
        if(n == 3) J3_on = 1.0;
    }

    // Decide whether to skip interactions. All the J2 values are taken later
    double              J1_rand = h5tb.param.J1_rand;
    std::vector<double> J2_rand = h5tb.param.J2_rand;
    double              J3_rand = h5tb.param.J3_rand;
    if(skip) {
        auto skip_this_pos = std::find(skip->begin(), skip->end(), get_position()) != skip->end();
        if(skip_this_pos) {
            // This site is skipped. No interactions should be contributed from here. Set all interactions to zero
            J1_rand = 0.0;
            J3_rand = 0.0;
            for(auto &J2r : J2_rand) J2r = 0;
        } else {
            // This site is not skipped, but we need to make sure not to interact with ones that are.
            for(const auto &r : J2_range) {
                if(r >= J2_rand.size()) break;
                auto skip_that_pos = std::find(skip->begin(), skip->end(), pos + r) != skip->end();
                if(skip_that_pos) J2_rand[r] = 0;
            }
            // Same for 3-body interactions
            auto skip_nnn1_pos = std::find(skip->begin(), skip->end(), pos + 1) != skip->end();
            auto skip_nnn2_pos = std::find(skip->begin(), skip->end(), pos + 2) != skip->end();
            if(skip_nnn1_pos or skip_nnn2_pos) J3_rand = 0;
        }
    }
    using namespace qm::spin::half;
    Eigen::Tensor<cplx, 4> MPO_nbody = MPO(); // Start with the full mpo
#pragma message "using sz instead of number operator"
    //    Eigen::Tensor<cplx, 2> n         = tenx::TensorCast(0.5 * (id + sz)); // Number operator
    Eigen::Tensor<cplx, 2> n = tenx::TensorMap(sz);
    Eigen::Tensor<cplx, 2> I = tenx::TensorMap(id); // identity

    for(const auto &r : J2_range) {
        double J2_count = 1.0;
        if(r >= J2_rand.size()) break;
        if(adjust_for_double_counting) {
            // Calculate double counting compensation
            // An interaction between sites i,j could be included multiple times in different multisite mpos.
            // Here we compensate for that so time evolution operators evolve the right amount
            // Example:
            //      Let L == 8 and R <= 3 and pos == 2
            //
            //      L               :  0,1,2,3,4,5,6,7
            //      multisite mpo[0]: [0,1,2,3]
            //      multisite mpo[1]:   [1,2,3,4]
            //      multisite mpo[2]:     [2,3,4,5]
            //      multisite mpo[3]:       [3,4,5,6]
            //      multisite mpo[4]:         [4,5,6,7]
            //
            //      Interaction counts for J[i,j], with i = posL and j = posL+r
            //      [posL,posL+1]   [posL,posL+2]  [posL,posL+3]
            //      [0,1]: 1        [0,2]: 1       [0,3]: 1
            //      [1,2]: 2        [1,3]: 2       [1,4]: 1
            //      [2,3]: 3        [2,4]: 2       [2,5]: 1
            //      [3,4]: 3        [3,5]: 2       [3,6]: 1
            //      [4,5]: 3        [4,6]: 2       [4,7]: 1
            //      [5,6]: 2        [5,7]: 1
            //      [6,7]: 1
            //
            //      Since in this example, pos == 2, we should compute the counts for [2,3], [2,4] and [2,5] in this function.
            //      If r == 2 in this loop, we should compute [2,4]: 2, because this mpo contributes J[2,4] twice as seen above.

            J2_count    = 0.0;     // For this particular r, for interaction from posL to posL+r
            size_t posI = pos;     // "i" in the interaction J(i,j)
            size_t posJ = pos + r; // "j" in the interaction J(i,j)
            for(size_t posL = 0; posL < settings::model::model_size - 1ul; posL++) {
                // posI and posJ are the left and right positions of a single 2-body interaction J(i,j).
                // posL and posR are the left and right edges of the multisite mpo
                size_t posR = posL + Rul;
                if(posR >= settings::model::model_size) break;
                if(posI >= posL and posJ <= posR) J2_count += 1.0; // Count if the interaction is in the multisite mpo
            }
            J2_count = std::max(J2_count, 1.0); // Avoid dividing by zero!
        }
        if(J2_count > 1.0) tools::log->trace("Adjusting for double counting: J2_count {} | pos {}", J2_count, get_position());
        J2_rand[r] /= J2_count;
    }
    MPO_nbody.slice(tenx::array4{F, 0, 0, 0}, extent4).reshape(extent2) = J1_on * (J1_rand * n - e_shift * I);

    if(R >= 1)
        for(const auto &r : J2_range) {
            if(r >= J2_rand.size()) break;
            MPO_nbody.slice(tenx::array4{F, static_cast<long>(r), 0, 0}, extent4).reshape(extent2) = J2_on * J2_rand[r] * n;
        }

    MPO_nbody.slice(tenx::array4{F, F - 1, 0, 0}, extent4).reshape(extent2) = J3_on * J3_rand * n;
    if constexpr(settings::debug_lbit)
        tools::log->info("Creating mpo | pos {} | nbody {} | skip {} \n{}", get_position(), nbody.value(), skip ? skip.value() : std::vector<size_t>{-1ul},
                         format_mpo(MPO_nbody, J1_rand, J2_rand, J3_rand));
    return MPO_nbody;
}

Eigen::Tensor<MpoSite::cplx, 4> LBit::MPO_shifted_view() const { return MPO_shifted_view(e_shift); }

Eigen::Tensor<MpoSite::cplx, 4> LBit::MPO_shifted_view(double site_energy) const {
    using namespace qm::spin::half;
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<cplx, 4> temp = MPO();
    long                   row  = temp.dimension(0) - 1;
    long                   col  = 0;
    if(parity_sep) row = temp.dimension(0) - 2;
#pragma message "using sz instead of number operator"

    //    Eigen::Tensor<cplx, 2> n                                           = tenx::TensorCast(0.5 * (id + sz));
    Eigen::Tensor<cplx, 2> n                                           = tenx::TensorMap(sz);
    Eigen::Tensor<cplx, 2> i                                           = tenx::TensorMap(id);
    temp.slice(tenx::array4{row, col, 0, 0}, extent4).reshape(extent2) = h5tb.param.J1_rand * n - site_energy * i;
    return temp;
}

std::unique_ptr<MpoSite> LBit::clone() const { return std::make_unique<LBit>(*this); }

long LBit::get_spin_dimension() const { return h5tb.param.spin_dim; }

void LBit::set_averages([[maybe_unused]] std::vector<TableMap> lattice_parameters, bool infinite) {
    tools::log->debug("LBIT MPO ({}): Setting averages", get_position());
    if(not infinite) {
        lattice_parameters.back()["J2_rand"]    = h5pp::varr_t<double>{0};
        lattice_parameters.back()["J3_rand"]    = 0.0;
        lattice_parameters.end()[-2]["J3_rand"] = 0.0;
    }
    if(parity_sep) throw std::runtime_error("Parity sector separation is not supported on LBIT MPO");
    set_parameters(lattice_parameters[get_position()]);
}

void LBit::save_hamiltonian(h5pp::File &file, std::string_view table_path) const {
    if(not file.linkExists(table_path)) file.createTable(h5tb.get_h5_type(), table_path, "LBIT");
    file.appendTableRecords(h5tb.param, table_path);
}

void LBit::load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) {
    auto ham_table = fmt::format("{}/hamiltonian", model_prefix);
    if(file.linkExists(ham_table)) {
        h5tb.param = file.readTableRecords<h5tb_lbit::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
    } else {
        throw except::runtime_error("Could not load MPO. Table [{}] does not exist", ham_table);
    }

    // Check that we are on the same point of the phase diagram
    using namespace settings::model::lbit;
    if(std::abs(h5tb.param.J1_mean - J1_mean) > 1e-6) throw except::runtime_error("J1_mean {:.16f} != {:.16f} lbit::J1_mean", h5tb.param.J1_mean, J1_mean);
    if(std::abs(h5tb.param.J2_mean - J2_mean) > 1e-6) throw except::runtime_error("J2_mean {:.16f} != {:.16f} lbit::J2_mean", h5tb.param.J2_mean, J2_mean);
    if(std::abs(h5tb.param.J3_mean - J3_mean) > 1e-6) throw except::runtime_error("J3_mean {:.16f} != {:.16f} lbit::J3_mean", h5tb.param.J3_mean, J3_mean);
    if(std::abs(h5tb.param.J1_wdth - J1_wdth) > 1e-6) throw except::runtime_error("J1_wdth {:.16f} != {:.16f} lbit::J1_wdth", h5tb.param.J1_wdth, J1_wdth);
    if(std::abs(h5tb.param.J2_wdth - J2_wdth) > 1e-6) throw except::runtime_error("J2_wdth {:.16f} != {:.16f} lbit::J2_wdth", h5tb.param.J2_wdth, J2_wdth);
    if(std::abs(h5tb.param.J3_wdth - J3_wdth) > 1e-6) throw except::runtime_error("J3_wdth {:.16f} != {:.16f} lbit::J3_wdth", h5tb.param.J3_wdth, J3_wdth);
    if(std::abs(h5tb.param.J2_xcls - J2_xcls) > 1e-6) throw except::runtime_error("J2_xcls {:.16f} != {:.16f} lbit::J2_xcls", h5tb.param.J2_xcls, J2_xcls);
    if(std::abs(h5tb.param.f_mixer - f_mixer) > 1e-6) throw except::runtime_error("f_mixer {:.16f} != {:.16f} lbit::f_mixer", h5tb.param.f_mixer, f_mixer);
    if(h5tb.param.distribution != distribution)
        throw except::runtime_error("distribution [{}] != lbit::distribution [{}]", h5tb.param.distribution, distribution);
    if(h5tb.param.J2_span != J2_span) throw except::runtime_error("J2_span {} != {} lbit::J2_span", h5tb.param.J2_span, J2_span);
    if(h5tb.param.u_layer != u_layer) throw except::runtime_error("u_layer {:.16f} != {:.16f} lbit::u_layer", h5tb.param.u_layer, u_layer);

    build_mpo();
}