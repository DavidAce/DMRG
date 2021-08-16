#include "EnvBase.h"
#include <config/debug.h>
#include <math/hash.h>
#include <math/num.h>
#include <math/tenx.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tools/common/log.h>
#include <utility>

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
EnvBase::EnvBase() : block(std::make_unique<Eigen::Tensor<Scalar, 3>>()) { assert_block(); } // default ctor
EnvBase::~EnvBase()                        = default;                                        // default dtor
EnvBase::EnvBase(EnvBase &&other) noexcept = default;                                        // default move ctor
EnvBase &EnvBase::operator=(EnvBase &&other) noexcept = default;                             // default move assign

EnvBase::EnvBase(const EnvBase &other)
    : edge_has_been_set(other.edge_has_been_set), block(std::make_unique<Eigen::Tensor<Scalar, 3>>(*other.block)), sites(other.sites), position(other.position),
      side(other.side), tag(other.tag), unique_id(other.unique_id), unique_id_mps(other.unique_id_mps), unique_id_mpo(other.unique_id_mpo),
      unique_id_env(other.unique_id_env) {
    assert_block();
}

EnvBase::EnvBase(size_t position_, std::string side_, std::string tag_)
    : block(std::make_unique<Eigen::Tensor<Scalar, 3>>()), position(position_), side(std::move(side_)), tag(std::move(tag_)) {
    assert_block();
}
EnvBase::EnvBase(std::string side_, std::string tag_, const MpsSite &MPS, const MpoSite &MPO)
    : block(std::make_unique<Eigen::Tensor<Scalar, 3>>()), side(std::move(side_)), tag(std::move(tag_)) {
    if(MPS.get_position() != MPO.get_position())
        throw std::logic_error(fmt::format("MPS and MPO have different positions: {} != {}", MPS.get_position(), MPO.get_position()));
    position = MPS.get_position();
    assert_block();
}

EnvBase &EnvBase::operator=(const EnvBase &other) {
    if(this != &other) {
        edge_has_been_set = other.edge_has_been_set;
        block             = std::make_unique<Eigen::Tensor<Scalar, 3>>(*other.block);
        sites             = other.sites;
        position          = other.position;
        side              = other.side;
        tag               = other.tag;
        unique_id         = other.unique_id;
        unique_id_mps     = other.unique_id_mps;
        unique_id_mpo     = other.unique_id_mpo;
        unique_id_env     = other.unique_id_env;
    }
    assert_block();
    return *this;
}

using Scalar = EnvBase::Scalar;

void EnvBase::assert_block() const {
    if(side != "L" and side != "R") throw std::runtime_error("Expected side [L|R]. Got: [" + side + "]");
    if(tag != "ene" and tag != "var") throw std::runtime_error("Expected tag [var|ene]. Got: [" + tag + "]");
    if(not block) throw std::runtime_error(fmt::format("env {} {} at pos {} is null", tag, side, get_position()));
}

void EnvBase::build_block(Eigen::Tensor<Scalar, 3> &otherblock, const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo) {
    /*!< Contracts a site into the block-> */
    // Note that otherblock, mps and mpo should correspond to the same site! I.e. their "get_position()" are all equal.
    // This can't be checked here though, so do that before calling this function.
    if constexpr(settings::debug) tools::log->trace("EnvBase::build_block(otherblock,mps,mpo): side({}), pos({})...", side, get_position());
    unique_id     = std::nullopt;
    unique_id_env = std::nullopt;
    unique_id_mps = std::nullopt;
    unique_id_mpo = std::nullopt;

    if(not block) block = std::make_unique<Eigen::Tensor<Scalar, 3>>();
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
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != mpo   dim[{}]:{}", side, tag, position.value(), 0,
                                                 mps.dimension(0), 2, mpo.dimension(2)));
        if(mps.dimension(1) != otherblock.dimension(0))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != left-block dim[{}]:{}", side, tag, position.value(), 1,
                                                 mps.dimension(1), 0, otherblock.dimension(0)));
        if(mpo.dimension(0) != otherblock.dimension(2))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mpo dim[{}]:{} != left-block dim[{}]:{}", side, tag, position.value(), 0,
                                                 mpo.dimension(0), 2, otherblock.dimension(2)));

        block->resize(mps.dimension(2), mps.dimension(2), mpo.dimension(1));
        block->device(tenx::omp::getDevice()) = otherblock.contract(mps, tenx::idx({0}, {1}))
                                                    .contract(mpo, tenx::idx({1, 2}, {0, 2}))
                                                    .contract(mps.conjugate(), tenx::idx({0, 3}, {1, 0}))
                                                    .shuffle(tenx::array3{0, 2, 1});
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
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != mpo dim[{}]:{}", side, tag, position.value(), 0,
                                                 mps.dimension(0), 2, mpo.dimension(2)));
        if(mps.dimension(2) != otherblock.dimension(0))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != right-block dim[{}]:{}", side, tag, position.value(), 2,
                                                 mps.dimension(2), 0, otherblock.dimension(0)));
        if(mpo.dimension(1) != otherblock.dimension(2))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mpo dim[{}]:{} != right-block dim[{}]:{}", side, tag, position.value(), 1,
                                                 mpo.dimension(1), 2, otherblock.dimension(2)));
        block->resize(mps.dimension(1), mps.dimension(1), mpo.dimension(0));
        block->device(tenx::omp::getDevice()) = otherblock.contract(mps, tenx::idx({0}, {2}))
                                                    .contract(mpo, tenx::idx({1, 2}, {1, 2}))
                                                    .contract(mps.conjugate(), tenx::idx({0, 3}, {2, 0}))
                                                    .shuffle(tenx::array3{0, 2, 1});
    }
    if constexpr(settings::debug) tools::log->trace("EnvBase::build_block(otherblock,mps,mpo): side({}), pos({})... OK", side, get_position());
}

void EnvBase::enlarge(const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo) {
    /*!< Contracts a site into the current block-> */

    // NOTE:
    // If side == "L", the mps and mpo should correspond to this env position+1
    // If side == "R", the mps and mpo should correspond to this env position-1
    // There is no way to check the positions here, so checks should be done before calling this function

    assert_block();
    Eigen::Tensor<Scalar, 3> thisblock = *block;
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
    // Do not clear the empty edge
    assert_block();
    if(sites > 0) {
        block->resize(0, 0, 0); // = Eigen::Tensor<Scalar,3>();
        tools::log->trace("Ejected env{} pos {}", side, get_position());
        unique_id = std::nullopt;
    }
    unique_id_env = std::nullopt;
    unique_id_mps = std::nullopt;
    unique_id_mpo = std::nullopt;
}

const Eigen::Tensor<Scalar, 3> &EnvBase::get_block() const {
    assert_block();
    return *block;
}
Eigen::Tensor<Scalar, 3> &EnvBase::get_block() {
    assert_block();
    return *block;
}

bool EnvBase::has_block() const { return block != nullptr and block->size() != 0; }

void EnvBase::assert_validity() const {
    assert_block();
    if(tenx::hasNaN(*block, fmt::format("env {} {}", tag, side))) {
        throw std::runtime_error(fmt::format("Environment {} side {} at position {} has NAN's", tag, side, get_position()));
    }
}

void EnvBase::assert_unique_id(const EnvBase &env, const MpsSite &mps, const MpoSite &mpo) const {
    std::string msg;
    if(env.get_unique_id() != unique_id_env) {
        msg.append(fmt::format("| env({}) {} !=", env.get_position(), env.get_unique_id()));
        if(unique_id_env) msg.append(fmt::format(" {} ", unique_id_env.value()));
    }
    if(mps.get_unique_id() != unique_id_mps) {
        msg.append(fmt::format("| mps({}) {} !=", mps.get_position(), mps.get_unique_id()));
        if(unique_id_mps) msg.append(fmt::format(" {} ", unique_id_mps.value()));
    }
    auto mpo_unique_id = tag == "ene" ? mpo.get_unique_id() : mpo.get_unique_id_sq();
    if(mpo_unique_id != unique_id_mpo) {
        msg.append(fmt::format("| mpo({}) {} !=", mpo.get_position(), mpo_unique_id));
        if(unique_id_mpo) msg.append(fmt::format(" {} ", unique_id_mpo.value()));
    }
    if(not msg.empty()) throw std::runtime_error(fmt::format("Environment {} side {}: unique id mismatch: {}", tag, side, msg));
}

bool EnvBase::is_real() const {
    assert_block();
    return tenx::isReal(*block, fmt::format("env {} {}", tag, side));
}

bool EnvBase::has_nan() const {
    assert_block();
    return tenx::hasNaN(*block, fmt::format("env {} {}", tag, side));
}

size_t EnvBase::get_position() const {
    if(position)
        return position.value();
    else
        throw std::runtime_error(fmt::format("Position hasn't been set on env side {}", side));
}

size_t EnvBase::get_sites() const { return sites; }

void EnvBase::set_edge_dims(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO, const Eigen::Tensor<Scalar, 1> &edge) {
    if(edge_has_been_set) return;
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
    sites             = 0;
    edge_has_been_set = true;
    unique_id         = std::nullopt;
    unique_id_env     = std::nullopt;
    unique_id_mps     = std::nullopt;
    unique_id_mpo     = std::nullopt;
}

std::size_t EnvBase::get_unique_id() const {
    if(unique_id) return unique_id.value();
    unique_id = hash::hash_buffer(get_block().data(), static_cast<size_t>(get_block().size()));
    return unique_id.value();
}

std::optional<std::size_t> EnvBase::get_unique_id_env() const { return unique_id_env; }
std::optional<std::size_t> EnvBase::get_unique_id_mps() const { return unique_id_mps; }
std::optional<std::size_t> EnvBase::get_unique_id_mpo() const { return unique_id_mpo; }

Eigen::Tensor<Scalar, 3> EnvBase::get_expansion_term(const MpsSite &mps, const MpoSite &mpo, double alpha) const {
    if constexpr(settings::debug)
        if(not num::all_equal(get_position(), mps.get_position(), mpo.get_position()))
            throw std::logic_error(fmt::format("class_env_{}::enlarge(): side({}), pos({}),: All positions are not equal: env {} | mps {} | mpo {}", tag, side,
                                               get_position(), get_position(), mps.get_position(), mpo.get_position()));

    Eigen::Tensor<Scalar, 4> mpo_tensor;
    if(tag == "ene")
        mpo_tensor = mpo.MPO();
    else if(tag == "var")
        mpo_tensor = mpo.MPO2();
    else
        throw std::runtime_error("Expected tag [var|ene]: Got: [" + tag + "]");

    if(side == "L") {
        long spin = mps.spin_dim();
        long chiL = mps.get_chiL();
        long chiR = mps.get_chiR() * mpo_tensor.dimension(1);
        if(get_block().dimension(0) != chiL) throw std::logic_error("block dim 0 != chiL");
        Eigen::Tensor<Scalar, 3> PL = get_block()
                                          .contract(mps.get_M_bare(), tenx::idx({0}, {1}))
                                          .contract(mpo_tensor, tenx::idx({1, 2}, {0, 2}))
                                          .shuffle(tenx::array4{3, 0, 1, 2})
                                          .reshape(tenx::array3{spin, chiL, chiR});
        PL = tenx::asNormalized(PL);
        return PL * PL.constant(alpha);

    } else if(side == "R") {
        long spin = mps.spin_dim();
        long chiL = mps.get_chiL() * mpo_tensor.dimension(0);
        long chiR = mps.get_chiR();
        if(get_block().dimension(0) != chiR) throw std::logic_error("block dim 0 != chiR");
        Eigen::Tensor<Scalar, 3> PR = get_block()
                                          .contract(mps.get_M_bare(), tenx::idx({0}, {2}))
                                          .contract(mpo_tensor, tenx::idx({1, 2}, {1, 2}))
                                          .shuffle(tenx::array4{3, 0, 1, 2})
                                          .reshape(tenx::array3{spin, chiL, chiR});
        PR = tenx::asNormalized(PR);
        return PR * PR.constant(alpha);

    } else
        throw std::runtime_error("Expected side L or R. Got: " + side);
}