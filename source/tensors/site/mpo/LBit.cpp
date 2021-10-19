#include "LBit.h"
#include <config/settings.h>
#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/linalg/tensor.h>
#include <math/num.h>
#include <math/rnd.h>
#include <math/tenx.h>
#include <qm/spin.h>

using namespace qm::spin::half;
using Scalar = std::complex<double>;

LBit::LBit(ModelType model_type_, size_t position_) : MpoSite(model_type_, position_) {
    h5tb.param.J1_mean  = settings::model::lbit::J1_mean;
    h5tb.param.J2_mean  = settings::model::lbit::J2_mean;
    h5tb.param.J3_mean  = settings::model::lbit::J3_mean;
    h5tb.param.J1_wdth  = settings::model::lbit::J1_wdth;
    h5tb.param.J2_wdth  = settings::model::lbit::J2_wdth;
    h5tb.param.J3_wdth  = settings::model::lbit::J3_wdth;
    h5tb.param.J2_xcls  = settings::model::lbit::J2_xcls;
    h5tb.param.J2_span  = settings::model::lbit::J2_span;
    h5tb.param.f_mixer  = settings::model::lbit::f_mixer;
    h5tb.param.u_layer  = settings::model::lbit::u_layer;
    h5tb.param.spin_dim = settings::model::lbit::spin_dim;

    // Adjust J2_span, it doesn't make sense to have it larger than the system size anyway, and setting it to 0
    // should automatically set it to system size
    if(h5tb.param.J2_span == 0 or h5tb.param.J2_span > settings::model::model_size) h5tb.param.J2_span = settings::model::model_size;

    copy_c_str(settings::model::lbit::distribution, h5tb.param.distribution);
    extent4 = {1, 1, h5tb.param.spin_dim, h5tb.param.spin_dim};
    extent2 = {h5tb.param.spin_dim, h5tb.param.spin_dim};
    h5tb_lbit::register_table_type();
    all_mpo_parameters_have_been_set =
        false; // There are no full lattice parameters on lbit, but we set it to false so we remember to call randomize on all sites in every model type
}

// double LBit::get_field() const { return h5tb.param.J1_rand; }
// double LBit::get_coupling() const { return h5tb.param.J2_rand + h5tb.param.J3_rand; }
void LBit::print_parameter_names() const { h5tb.print_parameter_names(); }
void LBit::print_parameter_values() const { h5tb.print_parameter_values(); }

void LBit::set_parameters(TableMap &parameters) {
    h5tb.param.J1_rand  = std::any_cast<double>(parameters["J1_rand"]);
    h5tb.param.J2_rand  = std::any_cast<h5tb_lbit::J2Type>(parameters["J2_rand"]);
    h5tb.param.J3_rand  = std::any_cast<double>(parameters["J3_rand"]);
    h5tb.param.J1_wdth  = std::any_cast<double>(parameters["J1_wdth"]);
    h5tb.param.J2_wdth  = std::any_cast<double>(parameters["J2_wdth"]);
    h5tb.param.J3_wdth  = std::any_cast<double>(parameters["J3_wdth"]);
    h5tb.param.J2_xcls  = std::any_cast<double>(parameters["J2_xcls"]);
    h5tb.param.f_mixer  = std::any_cast<double>(parameters["f_mixer"]);
    h5tb.param.u_layer  = std::any_cast<size_t>(parameters["u_layer"]);
    h5tb.param.spin_dim = std::any_cast<long>(parameters["spin_dim"]);
    std::strcpy(h5tb.param.distribution, std::any_cast<std::string>(parameters["distribution"]).c_str());
    all_mpo_parameters_have_been_set = true;
    build_mpo();
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
    parameters["f_mixer"]       = h5tb.param.f_mixer;
    parameters["u_layer"]       = h5tb.param.u_layer;
    parameters["spin_dim"]      = h5tb.param.spin_dim;
    parameters["distribution"]  = std::string(h5tb.param.distribution);
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
    if(not all_mpo_parameters_have_been_set) throw std::runtime_error("Improperly built MPO: Full lattice parameters haven't been set yet.");
    if(h5tb.param.J2_span >= h5tb.param.J2_rand.size())
        throw std::logic_error(fmt::format("expected J2_span ({}) < J2_rand.size()({})", h5tb.param.J2_span, h5tb.param.J2_rand.size()));
    tools::log->debug("LBIT MPO ({}): Building MPO", get_position());

    Eigen::Tensor<Scalar, 2> n = tenx::TensorCast(0.5 * (id + sz));
    Eigen::Tensor<Scalar, 2> I = tenx::TensorMap(id);
    long                     R = h5tb.param.J2_span;
    long                     F = R + 2l;
    mpo_internal.resize(F + 1, F + 1, h5tb.param.spin_dim, h5tb.param.spin_dim);
    mpo_internal.setZero();

    mpo_internal.slice(tenx::array4{0, 0, 0, 0}, extent4).reshape(extent2) = I;
    mpo_internal.slice(tenx::array4{F, F, 0, 0}, extent4).reshape(extent2) = I;
    mpo_internal.slice(tenx::array4{1, 0, 0, 0}, extent4).reshape(extent2) = n;

    if(R >= 2)
        for(const auto &i : num::range<long>(2, R + 1)) { mpo_internal.slice(tenx::array4{i, i - 1, 0, 0}, extent4).reshape(extent2) = I; }

    mpo_internal.slice(tenx::array4{F - 1, 1, 0, 0}, extent4).reshape(extent2) = n;
    mpo_internal.slice(tenx::array4{F, 0, 0, 0}, extent4).reshape(extent2)     = h5tb.param.J1_rand * n - e_reduced * I;

    if(R >= 1)
        for(const auto &i : num::range<long>(1, R + 1)) {
            mpo_internal.slice(tenx::array4{F, i, 0, 0}, extent4).reshape(extent2) = h5tb.param.J2_rand[static_cast<size_t>(i)] * n;
        }
    mpo_internal.slice(tenx::array4{F, F - 1, 0, 0}, extent4).reshape(extent2) = h5tb.param.J3_rand * n;

    if(tenx::hasNaN(mpo_internal)) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error(fmt::format("MPO at position {} has NAN's", get_position()));
    }
    unique_id = std::nullopt;
    build_mpo_squared();
}

Eigen::Tensor<Scalar, 1> LBit::get_MPO_edge_left() const {
    auto                     ldim = mpo_internal.dimension(0);
    Eigen::Tensor<Scalar, 1> ledge(ldim);
    ledge.setZero();
    ledge(ldim - 1) = 1;
    return ledge;
}

Eigen::Tensor<Scalar, 1> LBit::get_MPO_edge_right() const {
    auto                     rdim = mpo_internal.dimension(1);
    Eigen::Tensor<Scalar, 1> redge(rdim);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

Eigen::Tensor<Scalar, 1> LBit::get_MPO2_edge_left() const {
    auto edge = get_MPO_edge_left();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, tenx::idx()).reshape(tenx::array1{dim * dim});
}

Eigen::Tensor<Scalar, 1> LBit::get_MPO2_edge_right() const {
    auto edge = get_MPO_edge_right();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, tenx::idx()).reshape(tenx::array1{dim * dim});
}

void LBit::randomize_hamiltonian() {
    // J2(i,j) = exp(-|i-j|/J2_xcls) * Random_ij(J2_mean-J2_wdth_i/2, J2_mean+J2_wdth_i/2)
    using namespace settings::model::lbit;
    for(const auto &[r, J2_rand_ref] : iter::enumerate(h5tb.param.J2_rand)) {
        if(r > h5tb.param.J2_span) break;
        J2_rand_ref = std::exp(-static_cast<double>(r)/settings::model::lbit::J2_xcls); // exp(-r/J2_xcls) * Random_ij, where r = |i-j|, and Random_ij is compute d below
    }
    if(std::string(h5tb.param.distribution) == "normal") {
        h5tb.param.J1_rand = rnd::normal(J1_mean, J1_wdth);
        h5tb.param.J3_rand = rnd::normal(J3_mean, J3_wdth);
        for(auto &J2_rand_ref : h5tb.param.J2_rand) J2_rand_ref *= rnd::normal(J2_mean, J2_wdth);
    } else if(std::string(h5tb.param.distribution) == "lognormal") {
        h5tb.param.J1_rand = rnd::log_normal(J1_mean, J1_wdth);
        h5tb.param.J3_rand = rnd::log_normal(J3_mean, J3_wdth);
        for(auto &J2_rand_ref : h5tb.param.J2_rand) J2_rand_ref *= rnd::log_normal(J2_mean, J2_wdth);
    } else if(std::string(h5tb.param.distribution) == "uniform") {
        h5tb.param.J1_rand = rnd::uniform_double_box(J1_mean - J1_wdth / 2.0, J1_mean + J1_wdth / 2.0);
        h5tb.param.J3_rand = rnd::uniform_double_box(J3_mean - J3_wdth / 2.0, J3_mean + J3_wdth / 2.0);
        for(auto &J2_rand_ref : h5tb.param.J2_rand) J2_rand_ref *= rnd::uniform_double_box(J2_mean - J2_wdth / 2.0, J2_mean + J2_wdth / 2.0);
    } else if(std::string(h5tb.param.distribution) == "constant") {
        h5tb.param.J1_rand = settings::model::lbit::J1_mean;
        h5tb.param.J3_rand = settings::model::lbit::J3_mean;
        for(auto &J2_rand_ref : h5tb.param.J2_rand) J2_rand_ref *= settings::model::lbit::J2_mean;
    } else {
        throw std::runtime_error("Wrong distribution given. Expected one of <normal>, <lognormal>, <uniform> or <constant>");
    }
    all_mpo_parameters_have_been_set = false;
    mpo_squared                      = std::nullopt;
}

void LBit::set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) {
    switch(ptbMode) {
        case PerturbMode::ABSOLUTE: {
            h5tb.param.J1_pert = field_ptb;
            h5tb.param.J2_pert = coupling_ptb;
            h5tb.param.J3_pert = coupling_ptb;
            break;
        }
        case PerturbMode::PERCENTAGE: {
            h5tb.param.J1_pert = h5tb.param.J1_pert * field_ptb;
            h5tb.param.J2_pert = h5tb.param.J2_pert * coupling_ptb;
            h5tb.param.J3_pert = h5tb.param.J3_pert * coupling_ptb;
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_ABSOLUTE: {
            h5tb.param.J1_pert = rnd::uniform_double_box(-field_ptb, field_ptb);
            h5tb.param.J2_pert = rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.J3_pert = rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            break;
        }
        case PerturbMode::UNIFORM_RANDOM_PERCENTAGE: {
            h5tb.param.J1_pert = h5tb.param.J1_pert * rnd::uniform_double_box(-field_ptb, field_ptb);
            h5tb.param.J2_pert = h5tb.param.J2_pert * rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            h5tb.param.J3_pert = h5tb.param.J3_pert * rnd::uniform_double_box(-coupling_ptb, coupling_ptb);
            break;
        }
    }
    if(all_mpo_parameters_have_been_set) {
        Eigen::Tensor<Scalar, 2> n                                             = tenx::TensorCast(0.5 * (id + sz));
        Eigen::Tensor<Scalar, 2> I                                             = tenx::TensorMap(id);
        long                     F                                             = mpo_internal.dimension(0) - 1;
        mpo_internal.slice(tenx::array4{F, 0, 0, 0}, extent4).reshape(extent2) = h5tb.param.J1_rand * n - e_reduced * I;
        mpo_squared                                                            = std::nullopt;
        unique_id                                                              = std::nullopt;
        unique_id_sq                                                           = std::nullopt;
    }
    if(coupling_ptb == 0.0 and field_ptb == 0 and is_perturbed())
        throw std::runtime_error(fmt::format("MPO({}): Should have become unperturbed!", get_position()));
}

bool LBit::is_perturbed() const { return h5tb.param.J1_pert != 0.0 or h5tb.param.J2_pert != 0.0 or h5tb.param.J3_pert != 0.0; }

Eigen::Tensor<Scalar, 4> LBit::MPO_nbody_view(std::optional<std::vector<size_t>> nbody, std::optional<std::vector<size_t>> skip) const {
    // This function returns a view of the MPO including only n-body terms.
    // For instance, if nbody_terms == {2,3}, this would include 2-body and 3-body but exclude 1-body on-site terms.

    // Meanwhile, skip is a list of mpo positions that we should ignore. This lets us disable interactions with certain sites.
    // For instance if skip == {2}, then interaction terms such as J[0,2], J[1,2] and J[2,3] are set to zero.

    if(not nbody) return MPO();
    long   R                          = static_cast<long>(h5tb.param.J2_span); // Range
    long   F                          = R + 2l;                                // Final index of mpo
    size_t pos                        = get_position();
    auto   J2_range                   = num::range<size_t>(1, R + 1);
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
    double J1_rand = h5tb.param.J1_rand;
    double J3_rand = h5tb.param.J3_rand;
    auto   J2_rand = h5tb.param.J2_rand;
    if(skip) {
        for(const auto &s : skip.value()) {
            if(pos == s) {
                // This site is skipped. No interactions should be contributed from here
                J1_rand = 0.0;
                J3_rand = 0.0;
                for(auto &J2r : J2_rand) J2r = 0;
            } else {
                // This site is not skipped, but we need to make sure not to interact with the one that is skipped.
                for(const auto &r : J2_range) {
                    if(pos + r == s) J2_rand[r] = 0;
                }

                if(s == std::clamp<size_t>(s, pos, pos + 3)) J3_rand = 0;
            }
        }
    }

    Eigen::Tensor<Scalar, 4> MPO_nbody = MPO();                             // Start with the full mpo
    Eigen::Tensor<Scalar, 2> n         = tenx::TensorCast(0.5 * (id + sz)); // Number operator
    Eigen::Tensor<Scalar, 2> I         = tenx::TensorMap(id);               // identity

    auto J2_count = J2_rand;
    for(const auto &r : J2_range) {
        J2_count[r] = 1.0;
        if(adjust_for_double_counting) {
            // Calculate double counting compensation
            // An interaction between sites i,j could be included multiple times in different multisite mpos.
            // Here we compensate for that so time evolution operators evolve the right amount
            // Example:
            //      Let L == 8 and J2_span <= 3 and pos == 2
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

            J2_count[r] = 0.0;     // For this particular r, for interaction from posL to posL+r
            size_t posI = pos;     // "i" in the interaction J(i,j)
            size_t posJ = pos + r; // "j" in the interaction J(i,j)
            for(size_t posL = 0; posL < settings::model::model_size - 1ul; posL++) {
                // posI and posJ are the left and right positions of a single 2-body interaction J(i,j).
                // posL and posR are the left and right edges of the multisite mpo
                size_t posR = posL + h5tb.param.J2_span;
                if(posR >= settings::model::model_size) break;
                if(posI >= posL and posJ <= posR) J2_count[r] += 1.0; // Count if the interaction is in the multisite mpo
            }
            J2_count[r] = std::max(J2_count[r], 1.0); // Avoid dividing by zero!
        }
        J2_rand[r] /= J2_count[r];
    }

    MPO_nbody.slice(tenx::array4{F, 0, 0, 0}, extent4).reshape(extent2) = J1_on * (J1_rand * n - e_reduced * I);

    if(R >= 1)
        for(const auto &r : J2_range) { MPO_nbody.slice(tenx::array4{F, static_cast<long>(r), 0, 0}, extent4).reshape(extent2) = J2_on * J2_rand[r] * n; }

    MPO_nbody.slice(tenx::array4{F, F - 1, 0, 0}, extent4).reshape(extent2) = J3_on * J3_rand * n;

    return MPO_nbody;
}

Eigen::Tensor<Scalar, 4> LBit::MPO_reduced_view() const {
    if(e_reduced == 0) { return MPO(); }
    return MPO_reduced_view(e_reduced);
}

Eigen::Tensor<Scalar, 4> LBit::MPO_reduced_view(double site_energy) const {
    if(site_energy == 0) { return MPO(); }
    Eigen::Tensor<Scalar, 4> temp                                      = MPO();
    long                     row                                       = temp.dimension(0) - 1;
    long                     col                                       = 0;
    Eigen::Tensor<Scalar, 2> n                                         = tenx::TensorCast(0.5 * (id + sz));
    Eigen::Tensor<Scalar, 2> i                                         = tenx::TensorMap(id);
    temp.slice(tenx::array4{row, col, 0, 0}, extent4).reshape(extent2) = h5tb.param.J1_rand * n - site_energy * i;
    return temp;
}

std::unique_ptr<MpoSite> LBit::clone() const { return std::make_unique<LBit>(*this); }

long LBit::get_spin_dimension() const { return h5tb.param.spin_dim; }

void LBit::set_averages([[maybe_unused]] std::vector<TableMap> lattice_parameters, bool infinite, bool reverse) {
    tools::log->debug("LBIT MPO ({}): Setting averages", get_position());
    if(reverse) {
        std::reverse(lattice_parameters.begin(), lattice_parameters.end());
        for(size_t pos = 0; pos < lattice_parameters.size(); pos++) lattice_parameters[pos]["position"] = pos;
    }
    if(not infinite) {
        lattice_parameters.back()["J2_rand"]    = h5tb_lbit::J2Type{0};
        lattice_parameters.back()["J3_rand"]    = 0.0;
        lattice_parameters.end()[-2]["J3_rand"] = 0.0;
    }

    double J_sum = 0;
    for(auto &site_params : lattice_parameters) {
        auto J1_ = std::any_cast<double>(site_params["J1_rand"]);
        auto J3_ = std::any_cast<double>(site_params["J3_rand"]);
        auto J2_ = std::any_cast<h5tb_lbit::J2Type>(site_params["J2_rand"]);
        J_sum += J1_ + J3_;
        for(const auto &j2 : J2_) J_sum += j2;
    }
    if(parity_sep) psfactor = J_sum;
    set_parameters(lattice_parameters[get_position()]);
}

void LBit::save_hamiltonian(h5pp::File &file, std::string_view table_path) const {
    if(not file.linkExists(table_path)) file.createTable(h5tb_lbit::h5_type, table_path, "LBIT");
    file.appendTableRecords(h5tb, table_path);
    // Position 0 is also responsible for writing attributes
    if(position.value() != 0) return;

    file.writeAttribute(h5tb.param.J1_mean, "J1_mean", table_path);
    file.writeAttribute(h5tb.param.J2_mean, "J2_mean", table_path);
    file.writeAttribute(h5tb.param.J3_mean, "J3_mean", table_path);
    file.writeAttribute(h5tb.param.J1_wdth, "J1_wdth", table_path);
    file.writeAttribute(h5tb.param.J2_wdth, "J2_wdth", table_path);
    file.writeAttribute(h5tb.param.J3_wdth, "J3_wdth", table_path);
    file.writeAttribute(h5tb.param.J2_xcls, "J2_xcls", table_path);
    file.writeAttribute(h5tb.param.J2_span, "J2_span", table_path);
    file.writeAttribute(h5tb.param.f_mixer, "f_mixer", table_path);
    file.writeAttribute(h5tb.param.u_layer, "u_layer", table_path);
    file.writeAttribute(h5tb.param.distribution, "distribution", table_path);
    file.writeAttribute(h5tb.param.spin_dim, "spin_dim", table_path);
}

void LBit::load_hamiltonian(const h5pp::File &file, std::string_view model_prefix) {
    auto ham_table = fmt::format("{}/hamiltonian", model_prefix);
    if(file.linkExists(ham_table)) {
        h5tb.param                       = file.readTableRecords<h5tb_lbit::table>(ham_table, position);
        all_mpo_parameters_have_been_set = true;
        build_mpo();
    } else {
        throw std::runtime_error(fmt::format("Could not load MPO. Table [{}] does not exist", ham_table));
    }

    // Check that we are on the same point of the phase diagram
    using namespace settings::model::lbit;
    if(std::abs(h5tb.param.J1_mean - J1_mean) > 1e-6)
        throw std::runtime_error(fmt::format("J1_mean {:.16f} != {:.16f} lbit::J1_mean", h5tb.param.J1_mean, J1_mean));
    if(std::abs(h5tb.param.J2_mean - J2_mean) > 1e-6)
        throw std::runtime_error(fmt::format("J2_mean {:.16f} != {:.16f} lbit::J2_mean", h5tb.param.J2_mean, J2_mean));
    if(std::abs(h5tb.param.J3_mean - J3_mean) > 1e-6)
        throw std::runtime_error(fmt::format("J3_mean {:.16f} != {:.16f} lbit::J3_mean", h5tb.param.J3_mean, J3_mean));
    if(std::abs(h5tb.param.J1_wdth - J1_wdth) > 1e-6)
        throw std::runtime_error(fmt::format("J1_wdth {:.16f} != {:.16f} lbit::J1_wdth", h5tb.param.J1_wdth, J1_wdth));
    if(std::abs(h5tb.param.J2_wdth - J2_wdth) > 1e-6)
        throw std::runtime_error(fmt::format("J2_wdth {:.16f} != {:.16f} lbit::J2_wdth", h5tb.param.J2_wdth, J2_wdth));
    if(std::abs(h5tb.param.J3_wdth - J3_wdth) > 1e-6)
        throw std::runtime_error(fmt::format("J3_wdth {:.16f} != {:.16f} lbit::J3_wdth", h5tb.param.J3_wdth, J3_wdth));
    if(std::abs(h5tb.param.J2_xcls - J2_xcls) > 1e-6)
        throw std::runtime_error(fmt::format("J2_xcls {:.16f} != {:.16f} lbit::J2_xcls", h5tb.param.J2_xcls, J2_xcls));
    if(h5tb.param.J2_span != J2_span) throw std::runtime_error(fmt::format("J2_span {} != {} lbit::J2_span", h5tb.param.J2_span, J2_span));
    if(std::abs(h5tb.param.f_mixer - f_mixer) > 1e-6)
        throw std::runtime_error(fmt::format("f_mixer {:.16f} != {:.16f} lbit::f_mixer", h5tb.param.f_mixer, f_mixer));

    if(h5tb.param.u_layer != u_layer) throw std::runtime_error(fmt::format("u_layer {:.16f} != {:.16f} lbit::u_layer", h5tb.param.u_layer, u_layer));

    // We can use the mpo's on file here to check everything is correct
    std::string mpo_path = fmt::format("{}/mpo/H_{}", model_prefix, get_position());
    if(file.linkExists(mpo_path)) {
        auto mpo_dset = file.readDataset<Eigen::Tensor<Scalar, 4>>(mpo_path);
        if(tenx::VectorCast(MPO()) != tenx::VectorCast(mpo_dset)) throw std::runtime_error("Built MPO does not match the MPO on file");
    }
}