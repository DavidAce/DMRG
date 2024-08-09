#include "EnvBase.h"
#include "config/debug.h"
#include "debug/exceptions.h"
#include "math/hash.h"
#include "math/num.h"
#include "math/svd.h"
#include "math/tenx.h"
#include "math/tenx/threads.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tools/common/log.h"
#include <utility>
// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise, unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
EnvBase:: EnvBase() : block(std::make_unique<Eigen::Tensor<cplx, 3>>()) { assert_block(); } // default ctor
EnvBase::~EnvBase()                                    = default;                           // default dtor
EnvBase:: EnvBase(EnvBase &&other) noexcept            = default;                           // default move ctor
EnvBase  &EnvBase::operator=(EnvBase &&other) noexcept = default;                           // default move assign

EnvBase::EnvBase(const EnvBase &other)
    : block(std::make_unique<Eigen::Tensor<cplx, 3>>(*other.block)), sites(other.sites), position(other.position), side(other.side), tag(other.tag),
      unique_id(other.unique_id), unique_id_mps(other.unique_id_mps), unique_id_mpo(other.unique_id_mpo), unique_id_env(other.unique_id_env) {
    assert_block();
}

EnvBase::EnvBase(size_t position_, std::string side_, std::string tag_)
    : block(std::make_unique<Eigen::Tensor<cplx, 3>>()), position(position_), side(std::move(side_)), tag(std::move(tag_)) {
    assert_block();
}
EnvBase::EnvBase(std::string side_, std::string tag_, const MpsSite &MPS, const MpoSite &MPO)
    : block(std::make_unique<Eigen::Tensor<cplx, 3>>()), side(std::move(side_)), tag(std::move(tag_)) {
    if(MPS.get_position() != MPO.get_position())
        throw except::logic_error("MPS and MPO have different positions: {} != {}", MPS.get_position(), MPO.get_position());
    position = MPS.get_position();
    assert_block();
}

EnvBase &EnvBase::operator=(const EnvBase &other) {
    if(this != &other) {
        block         = std::make_unique<Eigen::Tensor<cplx, 3>>(*other.block);
        sites         = other.sites;
        position      = other.position;
        side          = other.side;
        tag           = other.tag;
        unique_id     = other.unique_id;
        unique_id_mps = other.unique_id_mps;
        unique_id_mpo = other.unique_id_mpo;
        unique_id_env = other.unique_id_env;
    }
    assert_block();
    return *this;
}

void EnvBase::assert_block() const {
    if(side != "L" and side != "R") throw except::runtime_error("Expected side [L|R]. Got: [{}]", side);
    if(tag != "ene" and tag != "var") throw except::runtime_error("Expected tag [var|ene]. Got: [{}]", tag);
    if(not block) throw except::runtime_error("env block {} {} at pos {} is nullptr", tag, side, get_position());
}

void EnvBase::build_block(Eigen::Tensor<cplx, 3> &otherblock, const Eigen::Tensor<cplx, 3> &mps, const Eigen::Tensor<cplx, 4> &mpo) {
    /*!< Contracts a site into the block-> */
    // Note that otherblock, mps and mpo should correspond to the same site! I.e. their "get_position()" are all equal.
    // This can't be checked here though, so do that before calling this function.
    if constexpr(settings::debug) tools::log->trace("EnvBase::build_block(otherblock,mps,mpo): side({}), pos({})", side, get_position());
    unique_id     = std::nullopt;
    unique_id_env = std::nullopt;
    unique_id_mps = std::nullopt;
    unique_id_mpo = std::nullopt;
    auto &threads = tenx::threads::get();
    if(not block) block = std::make_unique<Eigen::Tensor<cplx, 3>>();
    if(side == "L") {
        /*! # Left environment block contraction
         *   [       ]--0         [       ]--0 0--[LB]--1 1--[  GA    ]--2
         *   [       ]            [       ]                      |
         *   [       ]            [       ]                      0
         *   [       ]            [       ]
         *   [       ]            [       ]                      2
         *   [       ]            [       ]                      |
         *   [ block ]--2     =   [lblock ]--2            0--[   M    ]--1
         *   [       ]            [       ]                      |
         *   [       ]            [       ]                      3
         *   [       ]            [       ]
         *   [       ]            [       ]                      0
         *   [       ]            [       ]                      |
         *   [       ]--1         [       ]--1 0--[LB]--1  1--[GA conj ]--2
         */

        if(mps.dimension(0) != mpo.dimension(2))
            throw except::runtime_error("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != mpo   dim[{}]:{}", side, tag, position.value(), 0,
                                        mps.dimension(0), 2, mpo.dimension(2));
        if(mps.dimension(1) != otherblock.dimension(0))
            throw except::runtime_error("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != left-block dim[{}]:{}", side, tag, position.value(), 1,
                                        mps.dimension(1), 0, otherblock.dimension(0));
        if(mpo.dimension(0) != otherblock.dimension(2))
            throw except::runtime_error("env{} {} pos {} dimension mismatch: mpo dim[{}]:{} != left-block dim[{}]:{}", side, tag, position.value(), 0,
                                        mpo.dimension(0), 2, otherblock.dimension(2));
        // block->resize(mps.dimension(2), mps.dimension(2), mpo.dimension(1));
        // block->device(*threads->dev) = otherblock.contract(mps, tenx::idx({0}, {1}))
        //                                    .contract(mpo, tenx::idx({1, 2}, {0, 2}))
        //                                    .contract(mps.conjugate(), tenx::idx({0, 3}, {1, 0}))
        //                                    .shuffle(tenx::array3{0, 2, 1});
        Eigen::Tensor<cplx, 4> block_mps_mpo(otherblock.dimension(0), mps.dimension(2), mpo.dimension(1), mpo.dimension(2));
        block_mps_mpo.device(*threads->dev) = otherblock.contract(mps.conjugate(), tenx::idx({1}, {1})).contract(mpo, tenx::idx({1, 2}, {0, 3}));
        block->resize(mps.dimension(2), mps.dimension(2), mpo.dimension(1));
        block->device(*threads->dev) = mps.contract(block_mps_mpo, tenx::idx({0, 1}, {3, 0}));

    } else if(side == "R") {
        /*! # Right environment contraction
         *   0--[       ]          1--[   GB   ]--2  0--[LB]--1  0--[       ]
         *      [       ]                  |                        [       ]
         *      [       ]                  0                        [       ]
         *      [       ]                                           [       ]
         *      [       ]                  2                        [       ]
         *      [       ]                  |                        [       ]
         *   2--[ block ]     =     0--[   M   ]--1              2--[ rblock]
         *      [       ]                  |                        [       ]
         *      [       ]                  3                        [       ]
         *      [       ]                                           [       ]
         *      [       ]                  0                        [       ]
         *      [       ]                  |                        [       ]
         *   1--[       ]           1--[GB conj]--2  0--[LB]--1  1--[       ]
         */

        if(mps.dimension(0) != mpo.dimension(2))
            throw except::runtime_error("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != mpo dim[{}]:{}", side, tag, position.value(), 0,
                                        mps.dimension(0), 2, mpo.dimension(2));
        if(mps.dimension(2) != otherblock.dimension(0))
            throw except::runtime_error("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != right-block dim[{}]:{}", side, tag, position.value(), 2,
                                        mps.dimension(2), 0, otherblock.dimension(0));
        if(mpo.dimension(1) != otherblock.dimension(2))
            throw except::runtime_error("env{} {} pos {} dimension mismatch: mpo dim[{}]:{} != right-block dim[{}]:{}", side, tag, position.value(), 1,
                                        mpo.dimension(1), 2, otherblock.dimension(2));
        // block->resize(mps.dimension(1), mps.dimension(1), mpo.dimension(0));
        // block->device(*threads->dev) = otherblock.contract(mps, tenx::idx({0}, {2}))
        // .contract(mpo, tenx::idx({1, 2}, {1, 2}))
        // .contract(mps.conjugate(), tenx::idx({0, 3}, {2, 0}))
        // .shuffle(tenx::array3{0, 2, 1});

        Eigen::Tensor<cplx, 4> block_mps_mpo(otherblock.dimension(0), mps.dimension(1), mpo.dimension(0), mpo.dimension(2));
        block_mps_mpo.device(*threads->dev) = otherblock.contract(mps.conjugate(), tenx::idx({1}, {2})).contract(mpo, tenx::idx({1, 2}, {1, 3}));
        block->resize(mps.dimension(1), mps.dimension(1), mpo.dimension(0));
        block->device(*threads->dev) = mps.contract(block_mps_mpo, tenx::idx({0, 2}, {3, 0}));
    }
}

void EnvBase::enlarge(const Eigen::Tensor<cplx, 3> &mps, const Eigen::Tensor<cplx, 4> &mpo) {
    /*!< Contracts a site into the current block-> */

    // NOTE:
    // If side == "L", the mps and mpo should correspond to this env position+1
    // If side == "R", the mps and mpo should correspond to this env position-1
    // There is no way to check the positions here, so checks should be done before calling this function

    assert_block();
    Eigen::Tensor<cplx, 3> thisblock = *block;
    build_block(thisblock, mps, mpo);
    sites++;
    if(position) {
        if(side == "L")
            position = position.value() + 1;
        else if(position.value() > 0)
            position = position.value() - 1;
    }
}

void EnvBase::clear() {
    assert_block();
    block->resize(0, 0, 0); // = Eigen::Tensor<cplx,3>();
    tools::log->trace("Ejected env{} pos {}", side, get_position());
    unique_id = std::nullopt;

    // Do not clear the empty edge
    //    if(sites > 0) {
    //        block->resize(0, 0, 0); // = Eigen::Tensor<cplx,3>();
    //        tools::log->trace("Ejected env{} pos {}", side, get_position());
    //        unique_id = std::nullopt;
    //    }

    unique_id_env = std::nullopt;
    unique_id_mps = std::nullopt;
    unique_id_mpo = std::nullopt;
}

const Eigen::Tensor<cplx, 3> &EnvBase::get_block() const {
    assert_block();
    return *block;
}
Eigen::Tensor<cplx, 3> &EnvBase::get_block() {
    assert_block();
    return *block;
}

bool EnvBase::has_block() const { return block != nullptr and block->size() != 0; }

std::array<long, 3> EnvBase::get_dims() const { return get_block().dimensions(); }

void EnvBase::assert_validity() const {
    assert_block();
    if(tenx::hasNaN(*block)) { throw except::runtime_error("Environment {} side {} at position {} has NAN's", tag, side, get_position()); }
}

void EnvBase::assert_unique_id(const EnvBase &env, const MpsSite &mps, const MpoSite &mpo) const {
    std::vector<std::string> msg;
    if(unique_id_env and env.get_unique_id() != unique_id_env.value())
        msg.emplace_back(fmt::format("env({}): new {} | old {}", env.get_position(), env.get_unique_id(), unique_id_env.value()));
    if(unique_id_mps and mps.get_unique_id() != unique_id_mps.value())
        msg.emplace_back(fmt::format("mps({}): new {} | old {}", mps.get_position(), mps.get_unique_id(), unique_id_mps.value()));

    // mpo is special:
    // if this is an energy env we check the normal mpo, but
    // if this is a variance env we need to check the squared mpo
    if(tag == "ene") {
        if(unique_id_mpo and mpo.get_unique_id() != unique_id_mpo.value())
            msg.emplace_back(fmt::format("mpo({}): new {} | old {}", mpo.get_position(), mpo.get_unique_id(), unique_id_mpo.value()));
    } else if(tag == "var") {
        if(unique_id_mpo and mpo.get_unique_id_sq() != unique_id_mpo.value())
            msg.emplace_back(fmt::format("mpo_sq({}): new {} | old {}", mpo.get_position(), mpo.get_unique_id_sq(), unique_id_mpo.value()));
    } else {
        msg.emplace_back(fmt::format("Unrecognized tag: [{}]", tag));
    }
    if(not msg.empty())
        throw except::runtime_error("assert_unique_id: {}{}({}): unique id mismatch:\n{}\n"
                                    "Hint: remember to rebuild edges after operations that may modify them, "
                                    "like 1-site merge, move, normalization, projection, etc\n",
                                    tag, side, get_position(), fmt::join(msg, "\n"));
}

bool EnvBase::is_real() const {
    assert_block();
    return tenx::isReal(*block);
}

bool EnvBase::has_nan() const {
    assert_block();
    return tenx::hasNaN(*block);
}

size_t EnvBase::get_position() const {
    if(position)
        return position.value();
    else
        throw except::runtime_error("Position hasn't been set on env side {}", side);
}

size_t EnvBase::get_sites() const { return sites; }

void EnvBase::set_edge_dims(const Eigen::Tensor<cplx, 3> &MPS, const Eigen::Tensor<cplx, 4> &MPO, const Eigen::Tensor<cplx, 1> &edge) {
    assert_block();
    if(side == "L") {
        long mpsDim = MPS.dimension(1);
        long mpoDim = MPO.dimension(0);
        block->resize(tenx::array3{mpsDim, mpsDim, mpoDim});
        block->setZero();
        for(long i = 0; i < mpsDim; i++) {
            std::array<long, 1> extent1                     = {mpoDim};
            std::array<long, 3> offset3                     = {i, i, 0};
            std::array<long, 3> extent3                     = {1, 1, mpoDim};
            block->slice(offset3, extent3).reshape(extent1) = edge;
        }
    }
    if(side == "R") {
        long mpsDim = MPS.dimension(2);
        long mpoDim = MPO.dimension(1);
        block->resize(tenx::array3{mpsDim, mpsDim, mpoDim});
        block->setZero();
        for(long i = 0; i < mpsDim; i++) {
            std::array<long, 1> extent1                     = {mpoDim};
            std::array<long, 3> offset3                     = {i, i, 0};
            std::array<long, 3> extent3                     = {1, 1, mpoDim};
            block->slice(offset3, extent3).reshape(extent1) = edge;
        }
    }
    sites         = 0;
    unique_id     = std::nullopt;
    unique_id_env = std::nullopt;
    unique_id_mps = std::nullopt;
    unique_id_mpo = std::nullopt;
}

// void EnvBase::set_mixing_factor(double alpha) {
//     if(std::isnan(alpha)) throw except::logic_error("EnvBase::set_mixing_factor: alpha cannot be nan");
//     mixing_factor_alpha = alpha;
// }
// double EnvBase::get_mixing_factor() const { return mixing_factor_alpha; }

std::size_t EnvBase::get_unique_id() const {
    if(unique_id) return unique_id.value();
    unique_id = hash::hash_buffer(get_block().data(), safe_cast<size_t>(get_block().size()));
    return unique_id.value();
}

std::optional<std::size_t> EnvBase::get_unique_id_env() const { return unique_id_env; }
std::optional<std::size_t> EnvBase::get_unique_id_mps() const { return unique_id_mps; }
std::optional<std::size_t> EnvBase::get_unique_id_mpo() const { return unique_id_mpo; }

template<typename T>
Eigen::Tensor<T, 3> EnvBase::get_expansion_term(const MpsSite &mps, const MpoSite &mpo, double alpha, long rank_max) const {
    if constexpr(std::is_same_v<T, cplx>) {
        if(mps.is_real() and mpo.is_real()) { return get_expansion_term<real>(mps, mpo, alpha, rank_max).cast<cplx>(); }
    }

    assert(tag == "ene" or tag == "var");
    assert(side == "L" or side == "R");
    if constexpr(settings::debug)
        if(not num::all_equal(get_position(), mps.get_position(), mpo.get_position()))
            throw except::logic_error("class_env_{}::enlarge(): side({}), pos({}),: All positions are not equal: env {} | mps {} | mpo {}", tag, side,
                                      get_position(), get_position(), mps.get_position(), mpo.get_position());
    Eigen::Tensor<T, 4> mpo_tensor;
    Eigen::Tensor<T, 3> P;

    if constexpr(std::is_same_v<T, real>) {
        mpo_tensor = tag == "ene" ? mpo.MPO().real() : mpo.MPO2().real();
    } else {
        mpo_tensor = tag == "ene" ? mpo.MPO() : mpo.MPO2();
    }
    auto &threads = tenx::threads::get();

    if(side == "L") {
        // If P is too big we can pre-truncate here
        if(rank_max > 0 and mps.get_chiR() > std::max(rank_max, 32l)) {
            svd::solver svd;
            auto        mps_reduced = mps;
            // We can truncate the mps that goes into the expansion term here so that the SVD we do later doesn't become too large
            if constexpr(std::is_same_v<T, real>) {
                Eigen::Tensor<T, 3> M = mps_reduced.get_M().real();
                auto [U, S, V]        = svd.schmidt_into_left_normalized(M, mps_reduced.spin_dim(), svd::config(rank_max));
                if(mps_reduced.isCenter()) {
                    mps_reduced.set_M(U);
                    mps_reduced.set_LC(S);
                } else {
                    M.resize(U.dimensions());
                    M.device(*threads->dev) = U.contract(tenx::asDiagonal(S), tenx::idx({2}, {0}));
                    mps_reduced.set_M(M);
                }

            } else {
                auto [U, S, V] = svd.schmidt_into_left_normalized(mps_reduced.get_M(), mps_reduced.spin_dim(), svd::config(rank_max));
                if(mps_reduced.isCenter()) {
                    mps_reduced.set_M(U);
                    mps_reduced.set_LC(S);
                } else {
                    Eigen::Tensor<T, 3> M(U.dimensions());
                    M.device(*threads->dev) = U.contract(tenx::asDiagonal(S), tenx::idx({2}, {0}));
                    mps_reduced.set_M(M);
                }
            }
            tools::log->debug("Pre-truncated site {}: {} -> {} (rank max {})", mps.get_position(), mps.dimensions(), mps_reduced.dimensions(), rank_max);
            return get_expansion_term<T>(mps_reduced, mpo, alpha, -1);
        }

        long spin = mps.spin_dim();
        long chiL = mps.get_chiL();
        long chiR = mps.get_chiR() * mpo_tensor.dimension(1);
        assert(get_block().dimension(0) == chiL);
        P.resize(spin, chiL, chiR);
        if constexpr(std::is_same_v<T, real>) {
            Eigen::Tensor<T, 3> B   = get_block().real();
            Eigen::Tensor<T, 3> M   = mps.get_M().real();
            P.device(*threads->dev) = B.contract(M, tenx::idx({0}, {1}))
                                          .contract(mpo_tensor, tenx::idx({1, 2}, {0, 2}))
                                          .shuffle(tenx::array4{3, 0, 1, 2})
                                          .reshape(tenx::array3{spin, chiL, chiR});
        } else {
            P.device(*threads->dev) = get_block()
                                          .contract(mps.get_M(), tenx::idx({0}, {1}))
                                          .contract(mpo_tensor, tenx::idx({1, 2}, {0, 2}))
                                          .shuffle(tenx::array4{3, 0, 1, 2})
                                          .reshape(tenx::array3{spin, chiL, chiR});
        }

    } else if(side == "R") {
        // If P is too big we can pre-truncate here
        if(rank_max > 0 and mps.get_chiL() > std::max(rank_max, 32l)) {
            svd::solver svd;
            auto        mps_reduced = mps;
            // We can truncate the mps that goes into the expansion term here so that the SVD we do later doesn't become too large
            if constexpr(std::is_same_v<T, real>) {
                Eigen::Tensor<T, 3> M = mps_reduced.get_M().real();
                auto [U, S, V]        = svd.schmidt_into_right_normalized(M, mps_reduced.spin_dim(), svd::config(rank_max));
                M.resize(V.dimensions());
                M.device(*threads->dev) = tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(std::array{1, 0, 2});
                mps_reduced.set_M(M);
            } else {
                auto [U, S, V] = svd.schmidt_into_right_normalized(mps_reduced.get_M(), mps_reduced.spin_dim(), svd::config(rank_max));
                Eigen::Tensor<cplx, 3> M(V.dimensions());
                M.device(*threads->dev) = tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(std::array{1, 0, 2});
                mps_reduced.set_M(M);
            }
            tools::log->debug("Pre-truncated site {}: {} -> {} (rank max {})", mps.get_position(), mps.dimensions(), mps_reduced.dimensions(), rank_max);
            return get_expansion_term<T>(mps_reduced, mpo, alpha, -1);
        }

        long spin = mps.spin_dim();
        long chiL = mps.get_chiL() * mpo_tensor.dimension(0);
        long chiR = mps.get_chiR();
        assert(get_block().dimension(0) == chiR);
        P.resize(spin, chiL, chiR);
        if constexpr(std::is_same_v<T, real>) {
            Eigen::Tensor<T, 3> B   = get_block().real();
            Eigen::Tensor<T, 3> M   = mps.get_M().real();
            P.device(*threads->dev) = B.contract(M, tenx::idx({0}, {2}))
                                          .contract(mpo_tensor, tenx::idx({1, 2}, {1, 2}))
                                          .shuffle(tenx::array4{3, 0, 1, 2})
                                          .reshape(tenx::array3{spin, chiL, chiR});
        } else {
            P.device(*threads->dev) = get_block()
                                          .contract(mps.get_M(), tenx::idx({0}, {2}))
                                          .contract(mpo_tensor, tenx::idx({1, 2}, {1, 2}))
                                          .shuffle(tenx::array4{3, 0, 1, 2})
                                          .reshape(tenx::array3{spin, chiL, chiR});
        }
    }
    // P         = tenx::asNormalized(P);
    // if(alpha  <= 0.0) alpha = get_mixing_factor();
    if(std::isnan(alpha)) throw except::logic_error("the mixing factor for env {}{}({}) is nan", tag, side, get_position());
    if(alpha <= 0) throw except::logic_error("the mixing factor for env {}{}({}) is not positive: {:.3e}", tag, side, get_position(), alpha);
    return P * P.constant(alpha);
}

template Eigen::Tensor<real, 3> EnvBase::get_expansion_term<real>(const MpsSite &mps, const MpoSite &mpo, double alpha, long rank_max) const;
template Eigen::Tensor<cplx, 3> EnvBase::get_expansion_term<cplx>(const MpsSite &mps, const MpoSite &mpo, double alpha, long rank_max) const;
