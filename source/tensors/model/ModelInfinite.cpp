#include "ModelInfinite.h"
#include "tensors/site/mpo/MpoFactory.h"
#include <config/settings.h>
#include <general/iter.h>
#include <math/num.h>
#include <math/svd.h>
#include <math/tenx.h>

ModelInfinite::ModelInfinite() = default;

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
ModelInfinite::~ModelInfinite()                     = default;            // default dtor
ModelInfinite::ModelInfinite(ModelInfinite &&other) = default;            // default move ctor
ModelInfinite &ModelInfinite::operator=(ModelInfinite &&other) = default; // default move assign

ModelInfinite::ModelInfinite(const ModelInfinite &other) : cache(other.cache), HA(other.HA->clone()), HB(other.HB->clone()), model_type(other.model_type) {}

ModelInfinite &ModelInfinite::operator=(const ModelInfinite &other) {
    // check for self-assignment
    if(this != &other) {
        cache      = other.cache;
        HA         = other.HA->clone();
        HB         = other.HB->clone();
        model_type = other.model_type;
    }
    return *this;
}

void ModelInfinite::initialize(ModelType model_type_) {
    tools::log->trace("Initializing model");
    // Generate MPO
    model_type = model_type_;
    HA         = MpoFactory::create_mpo(0, model_type);
    HB         = MpoFactory::create_mpo(1, model_type);
}

void ModelInfinite::randomize() {
    tools::log->trace("Randomizing hamiltonian");
    HA->randomize_hamiltonian();
    HB->randomize_hamiltonian();
    std::vector<MpoSite::TableMap> all_params;
    all_params.push_back(HA->get_parameters());
    all_params.push_back(HB->get_parameters());
    HA->set_averages(all_params, true);
    HB->set_averages(all_params, true);
}

void ModelInfinite::reset_mpo_squared() {
    tools::log->debug("Resetting squared MPO");
    HA->build_mpo_squared();
    HB->build_mpo_squared();
}

void ModelInfinite::rebuild_mpo_squared(std::optional<SVDMode> svdMode) {
    if(settings::precision::use_compressed_mpo_squared_all) {
        tools::log->trace("Compressing MPO²");
        throw std::runtime_error("Compressing the squared MPO² is currently unsupported on infinite systems.\n"
                                 "Set settings::precision::use_compressed_mpo_squared_all = false");
        auto mpo_compressed = get_compressed_mpo_squared(svdMode);
        HA->set_mpo_squared(mpo_compressed[0]);
        HA->set_mpo_squared(mpo_compressed[1]);
    } else
        reset_mpo_squared();
}

std::vector<Eigen::Tensor<ModelInfinite::Scalar, 4>> ModelInfinite::get_compressed_mpo_squared(std::optional<SVDMode> svdMode) {
    // First, rebuild the MPO's
    std::vector<Eigen::Tensor<Scalar, 4>> mpos_sq;
    mpos_sq.emplace_back(HA->get_non_compressed_mpo_squared());
    mpos_sq.emplace_back(HB->get_non_compressed_mpo_squared());

    // Setup SVD
    // Here we need a lot of precision:
    //  - Use very low svd threshold
    //  - Force the use of JacobiSVD by setting the switchsize to something large
    //  - Force the use of Lapacke -- it is more precise than Eigen (I don't know why)
    svd::solver svd;
    svd.threshold   = 1e-12;
    svd.switchsize  = 50000;
    svd.use_lapacke = true;
    if(svdMode == SVDMode::EIGEN) svd.use_lapacke = false;

    // Print the results
    std::vector<std::string> report;
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &mpo : mpos_sq) report.emplace_back(fmt::format("{}", mpo.dimensions()));

    for(size_t iter = 0; iter < 1; iter++) {
        // Next compress from left to right
        Eigen::Tensor<Scalar, 2> T_l2r; // Transfer matrix
        Eigen::Tensor<Scalar, 4> T_mpo_sq;
        for(const auto &[idx, mpo_sq] : iter::enumerate(mpos_sq)) {
            if(T_l2r.size() == 0)
                T_mpo_sq = mpo_sq;
            else
                T_mpo_sq = T_l2r.contract(mpo_sq, tenx::idx({1}, {0}));

            if(idx == mpos_sq.size() - 1) {
                mpo_sq = T_mpo_sq;
            } else {
                auto [U, S, V] = svd.split_mpo_l2r(T_mpo_sq);
                T_l2r          = tenx::asDiagonal(S).contract(V, tenx::idx({1}, {0}));
                if(idx < mpos_sq.size() - 1)
                    mpo_sq = U;
                else
                    // The remaining transfer matrix T can be multiplied back into the last MPO from the right
                    mpo_sq = U.contract(T_l2r, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
            }
        }

        // Now we have done left to right. Next we do right to left
        Eigen::Tensor<Scalar, 2> T_r2l;    // Transfer matrix
        Eigen::Tensor<Scalar, 4> mpo_sq_T; // Absorbs transfer matrix
        for(const auto &[idx, mpo_sq] : iter::enumerate_reverse(mpos_sq)) {
            if(T_r2l.size() == 0)
                mpo_sq_T = mpo_sq;
            else
                mpo_sq_T = mpo_sq.contract(T_r2l, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
            if(idx == 0) {
                mpo_sq = mpo_sq_T;
            } else {
                auto [U, S, V] = svd.split_mpo_r2l(mpo_sq_T);
                T_r2l          = U.contract(tenx::asDiagonal(S), tenx::idx({1}, {0}));
                if(idx > 0)
                    mpo_sq = V;
                else
                    // The remaining transfer matrix T can be multiplied back into the first MPO from the left
                    mpo_sq = T_r2l.contract(V, tenx::idx({1}, {0}));
            }
        }
    }

    // Print the results
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &[idx, msg] : iter::enumerate(report)) tools::log->trace("MPO² {}: {} -> {}", idx, msg, mpos_sq[idx].dimensions());
    return mpos_sq;
}

bool ModelInfinite::is_real() const { return HA->is_real() and HB->is_real(); }
bool ModelInfinite::has_nan() const { return HA->has_nan() or HB->has_nan(); }
void ModelInfinite::assert_validity() const {
    HA->assert_validity();
    HB->assert_validity();
}

const MpoSite &ModelInfinite::get_mpo_siteA() const { return *HA; }
const MpoSite &ModelInfinite::get_mpo_siteB() const { return *HB; }
MpoSite       &ModelInfinite::get_mpo_siteA() { return *HA; }
MpoSite       &ModelInfinite::get_mpo_siteB() { return *HB; }

Eigen::DSizes<long, 4> ModelInfinite::dimensions() const {
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    return Eigen::DSizes<long, 4>{dim0, dim1, dim2, dim3};
}

bool ModelInfinite::is_reduced() const { return HA->is_reduced() and HB->is_reduced(); }

double ModelInfinite::get_energy_per_site_reduced() const {
    if(not num::all_equal(HA->get_reduced_energy(), HB->get_reduced_energy()))
        throw std::runtime_error(fmt::format("Reduced energy mismatch: HA {:.16f} != HB {:.16f}", HA->get_reduced_energy(), HB->get_reduced_energy()));
    return HA->get_reduced_energy();
}

void ModelInfinite::set_reduced_energy_per_site(double site_energy) {
    HA->set_reduced_energy(site_energy);
    HB->set_reduced_energy(site_energy);
}

const Eigen::Tensor<ModelInfinite::Scalar, 4> &ModelInfinite::get_2site_mpo_AB() const {
    if(cache.twosite_mpo_AB) return cache.twosite_mpo_AB.value();
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    cache.twosite_mpo_AB =
        HA->MPO().contract(HB->MPO(), tenx::idx({1}, {0})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(tenx::array4{dim0, dim1, dim2, dim3});
    return cache.twosite_mpo_AB.value();
}

const Eigen::Tensor<ModelInfinite::Scalar, 4> &ModelInfinite::get_2site_mpo_BA() const {
    if(cache.twosite_mpo_BA) return cache.twosite_mpo_BA.value();
    long dim0 = HB->MPO().dimension(0);
    long dim1 = HA->MPO().dimension(1);
    long dim2 = HB->MPO().dimension(2) * HA->MPO().dimension(2);
    long dim3 = HB->MPO().dimension(3) * HA->MPO().dimension(3);
    cache.twosite_mpo_BA =
        HB->MPO().contract(HA->MPO(), tenx::idx({1}, {0})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(tenx::array4{dim0, dim1, dim2, dim3});
    return cache.twosite_mpo_BA.value();
}

const Eigen::Tensor<ModelInfinite::Scalar, 2> &ModelInfinite::get_2site_ham_AB() const {
    if(cache.twosite_ham_AB) return cache.twosite_ham_AB.value();
    auto twosite_mpo_AB  = get_2site_mpo_AB();
    auto edgeL           = get_mpo_siteA().get_MPO_edge_left();
    auto edgeR           = get_mpo_siteB().get_MPO_edge_right();
    cache.twosite_ham_AB = twosite_mpo_AB.contract(edgeL, tenx::idx({0}, {0})).contract(edgeR, tenx::idx({0}, {0}));
    return cache.twosite_ham_AB.value();
}
const Eigen::Tensor<ModelInfinite::Scalar, 2> &ModelInfinite::get_2site_ham_BA() const {
    if(cache.twosite_ham_BA) return cache.twosite_ham_BA.value();
    auto twosite_mpo_BA  = get_2site_mpo_BA();
    auto edgeL           = get_mpo_siteB().get_MPO_edge_left();
    auto edgeR           = get_mpo_siteA().get_MPO_edge_right();
    cache.twosite_ham_BA = twosite_mpo_BA.contract(edgeL, tenx::idx({0}, {0})).contract(edgeR, tenx::idx({0}, {0}));
    return cache.twosite_ham_BA.value();
}

void ModelInfinite::clear_cache() { cache = Cache(); }
