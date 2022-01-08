

#include "qm/gate.h"
#include <Eigen/Core>
#include <general/iter.h>
#include <math/linalg/tensor.h>
#include <math/num.h>
#include <math/tenx.h>
#include <set>
#include <tools/common/log.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <utility>

template<typename T>
std::vector<T> subset(const std::vector<T> &vec, size_t idx_start, size_t num) {
    if(idx_start + num > vec.size())
        throw std::range_error(fmt::format("Vector subset start {} num {} out of range for vector of size {}", idx_start, num, vec.size()));
    auto vec_bgn = vec.begin() + static_cast<long>(idx_start);
    auto vec_end = vec_bgn + static_cast<long>(num);
    return std::vector<T>(vec_bgn, vec_end);
}

template<auto N, typename T, auto M>
std::array<T, N> subset(const std::array<T, M> &arr, size_t idx_start) {
    if(idx_start + N > arr.size())
        throw std::range_error(fmt::format("Vector subset start {} num {} out of range for vector of size {}", idx_start, N, arr.size()));
    auto             vec_bgn = arr.begin() + idx_start;
    auto             vec_end = vec_bgn + N;
    std::array<T, N> res{};
    std::copy(vec_bgn, vec_end, res.begin());
    return res;
}

template<typename T>
T prod(const std::vector<T> &vec, size_t idx_start, size_t num, T seed = 1) {
    auto sub = subset(vec, idx_start, num);
    return std::accumulate(sub.begin(), sub.end(), std::multiplies(seed));
}

template<auto N, typename T, auto M>
T prod(const std::array<T, M> &arr, size_t idx_start, T seed = 1) {
    auto sub = subset<N>(arr, idx_start);
    T    acc = seed;
    for(auto &s : sub) acc *= s;
    return acc;
}

template<typename T>
std::vector<T> concat(const std::vector<T> &v1, const std::vector<T> &v2) {
    auto res = v1;
    res.insert(res.end(), v2.begin(), v2.end());
    return res;
}

template<typename T, auto N, auto M>
auto concat(const std::array<T, N> &a1, const std::array<T, M> &a2) {
    std::array<T, N + M> result{};
    std::copy(a1.cbegin(), a1.cend(), result.begin());
    std::copy(a2.cbegin(), a2.cend(), result.begin() + N);
    return result;
}

template<typename T, auto N>
auto repeat(const std::array<T, N> &a) {
    return concat(a, a);
}

auto group(const std::vector<long> &dim, const std::vector<size_t> &pattern) {
    std::vector<long> res(pattern.size());
    size_t            dim_offset = 0;
    for(size_t i = 0; i < pattern.size(); i++) {
        long product = 1;
        for(size_t j = 0; j < pattern[i]; j++) product *= dim[dim_offset + j];
        res[i] = product;
        dim_offset += pattern[i];
    }
    return res;
}

template<auto N, auto M>
auto group(const std::array<long, N> &dim, const std::array<size_t, M> &pattern) {
    // Group tensor indices together.
    std::array<long, M> res{};
    size_t              dim_offset = 0;
    for(size_t i = 0; i < M; i++) {
        long product = 1;
        for(size_t j = 0; j < pattern[i]; j++) product *= dim[dim_offset + j];
        res[i] = product;
        dim_offset += pattern[i];
    }
    return res;
}

template<typename T1, typename T2>
void erase(std::vector<T1> &vec, T2 val) {
    auto it = vec.begin();
    while(it != vec.end()) {
        if(*it == static_cast<T1>(val))
            it = vec.erase(it);
        else
            ++it;
    }
}

Eigen::Tensor<qm::cplx, 2> contract_a(const Eigen::Tensor<qm::cplx, 2> &m, const Eigen::Tensor<qm::cplx, 2> &ud, const std::array<long, 4> &shp_mid4,
                                      const std::array<long, 4> &shp_udn4, const std::array<long, 6> &shf6, const tenx::idxlistpair<1> &idx1,
                                      const tenx::idxlistpair<2> &idx2, const std::array<long, 2> &dim2) {
    return ud.reshape(shp_udn4).contract(m.reshape(shp_mid4), idx1).contract(ud.conjugate().reshape(shp_udn4), idx2).shuffle(shf6).reshape(dim2);
}

Eigen::Tensor<qm::cplx, 2> contract_b(const Eigen::Tensor<qm::cplx, 2> &m, const Eigen::Tensor<qm::cplx, 2> &ud, const std::array<long, 2> &shp_udn2,
                                      const std::array<long, 4> &shp_udn4, const tenx::idxlistpair<1> &idx1, const tenx::idxlistpair<2> &idx2) {
    return ud.reshape(shp_udn4).contract(m, idx1).contract(ud.conjugate().reshape(shp_udn4), idx2).reshape(shp_udn2);
}

Eigen::Tensor<qm::cplx, 2> contract_c(const Eigen::Tensor<qm::cplx, 2> &m, const Eigen::Tensor<qm::cplx, 2> &ud, const std::array<long, 6> &shp_mid6,
                                      const tenx::idxlistpair<1> &idx_up, const tenx::idxlistpair<1> &idx_dn, const std::array<long, 6> &shf6,
                                      const std::array<long, 2> &dim2) {
    return ud.contract(m.reshape(shp_mid6), idx_up).contract(ud.conjugate(), idx_dn).shuffle(shf6).reshape(dim2);
}

Eigen::Tensor<qm::cplx, 2> contract_d(const Eigen::Tensor<qm::cplx, 2> &m, const Eigen::Tensor<qm::cplx, 2> &ud, const std::array<long, 4> &shp_mid4,
                                      const tenx::idxlistpair<1> &idx_up, const tenx::idxlistpair<1> &idx_dn, const std::array<long, 4> &shf4,
                                      const std::array<long, 2> &dim2) {
    return ud.contract(m.reshape(shp_mid4), idx_up).contract(ud.conjugate(), idx_dn).shuffle(shf4).reshape(dim2);
}

Eigen::Tensor<qm::cplx, 2> qm::Gate::exp_internal(const Eigen::Tensor<cplx, 2> &op_, cplx alpha) const {
    /* Note for flbit:
     *  Let h = op(i,i), i.e. h are the diagonal entries in op
     *  Let alpha = -i * delta be purely imaginary (since delta is real). So delta = imag(-alpha)
     *  Now notice that
     *       exp( -i * delta * h )
     *  can become imprecise when delta is large, since delta and h are real and h ~O(1)
     *
     *  Instead, in the flbit case, we can exploit that h is real to compute
     *       exp(-i * mod(delta * real(h), 2*pi))
     *  which is equivalent but with a much smaller exponent.
     *
     *  Remember to do the modulo separately on each diagonal entry h!
     */
    auto op_map = tenx::MatrixMap(op_);
    if(op_map.isDiagonal() and op_map.imag().isZero() and std::real(alpha) == 0) {
        auto minus_i = std::complex<double>(0, -1);

        auto diag =
            op_map.diagonal()
                .unaryViewExpr([&alpha, &minus_i](const cplx &h) { return std::exp(minus_i * std::fmod(std::imag(-alpha) * std::real(h), 2.0 * M_PI)); })
                .asDiagonal();
        return tenx::TensorMap(diag.toDenseMatrix());
    } else {
        return tenx::TensorMap((alpha * tenx::MatrixMap(op_)).exp().eval());
    }
}

qm::Gate::Gate(const Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic> &op_, std::vector<size_t> pos_, std::vector<long> dim_)
    : pos(std::move(pos_)), dim(std::move(dim_)) {
    auto dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1l, std::multiplies<>());
    if(dim_prod != op_.rows() or dim_prod != op_.cols())
        throw std::logic_error(fmt::format("dim {} not compatible with matrix dimensions {} x {}", dim, op_.rows(), op_.cols()));
    if(pos.size() != dim.size()) throw std::logic_error(fmt::format("pos.size() {} != dim.size() {}", pos, dim));
    op = tenx::TensorMap(op_);
}

qm::Gate::Gate(const Eigen::Tensor<cplx, 2> &op_, std::vector<size_t> pos_, std::vector<long> dim_) : op(op_), pos(std::move(pos_)), dim(std::move(dim_)) {
    auto dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies<>());
    if(dim_prod != op_.dimension(0) or dim_prod != op_.dimension(1))
        throw std::logic_error(fmt::format("dim {} not compatible with matrix dimensions {} x {}", dim, op_.dimension(0), op_.dimension(1)));
    if(pos.size() != dim.size()) throw std::logic_error(fmt::format("pos.size() {} != dim.size() {}", pos, dim));
}
qm::Gate::Gate(const Eigen::Tensor<cplx, 2> &op_, std::vector<size_t> pos_, std::vector<long> dim_, cplx alpha) : pos(std::move(pos_)), dim(std::move(dim_)) {
    auto dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies<>());
    if(dim_prod != op_.dimension(0) or dim_prod != op_.dimension(1))
        throw std::logic_error(fmt::format("dim {} not compatible with matrix dimensions {} x {}", dim, op_.dimension(0), op_.dimension(1)));
    if(pos.size() != dim.size()) throw std::logic_error(fmt::format("pos.size() {} != dim.size() {}", pos, dim));
    op = exp_internal(op_, alpha);
}

void                        qm::Gate::mark_as_used() const { used = true; }
void                        qm::Gate::unmark_as_used() const { used = false; }
bool                        qm::Gate::was_used() const { return used; }
void                        qm::Gate::exp_inplace(cplx alpha) { op = exp_internal(op, alpha); }
qm::Gate                    qm::Gate::exp(cplx alpha) const { return Gate(op, pos, dim, alpha); }
bool                        qm::Gate::isUnitary(double prec) const { return tenx::MatrixMap(op).isUnitary(prec); }
Eigen::Tensor<qm::cplx, 2> &qm::Gate::adjoint() const {
    if(adj) return adj.value();
    adj = op.conjugate().shuffle(tenx::array2{1, 0});
    return adj.value();
}

template<auto rank>
std::array<long, rank> qm::Gate::shape() const {
    if constexpr(rank == 1) { return {op.dimension(0) * op.dimension(1)}; }
    if constexpr(rank == 2) {
        return {op.dimension(0), op.dimension(1)};
    } else {
        if(rank == 2 * dim.size()) {
            std::array<long, rank> dims{};
            for(size_t i = 0; i < dims.size(); i++) dims[i] = dim[i % dim.size()];
            return dims;
        } else {
            throw std::range_error(fmt::format("Can't compute shape of rank {} for gate with pos {} and dim {}", rank, pos, dim));
        }
    }
}
template std::array<long, 1>  qm::Gate::shape() const;
template std::array<long, 2>  qm::Gate::shape() const;
template std::array<long, 4>  qm::Gate::shape() const;
template std::array<long, 6>  qm::Gate::shape() const;
template std::array<long, 8>  qm::Gate::shape() const;
template std::array<long, 10> qm::Gate::shape() const;
template std::array<long, 12> qm::Gate::shape() const;
template std::array<long, 14> qm::Gate::shape() const;
template std::array<long, 16> qm::Gate::shape() const;

std::vector<size_t> qm::Gate::idx() const {
    std::vector<size_t> idx(pos.size() * 2);
    std::iota(idx.begin(), idx.end(), 0);
    return idx;
}

std::vector<size_t> qm::Gate::idx_up() const {
    std::vector<size_t> idx(pos.size());
    std::iota(idx.begin(), idx.end(), 0);
    return idx;
}

std::vector<size_t> qm::Gate::idx_dn() const {
    std::vector<size_t> idx(pos.size());
    std::iota(idx.begin(), idx.end(), pos.size());
    return idx;
}

std::vector<size_t> qm::Gate::idx(const std::vector<size_t> &pos_) const {
    std::vector<size_t> idx;
    for(const auto &p_ : pos_)
        for(const auto &[i, p] : iter::enumerate(pos))
            if(p == p_) {
                idx.emplace_back(i);
                idx.emplace_back(i + pos.size());
            }
    return idx;
}

std::vector<size_t> qm::Gate::idx_up(const std::vector<size_t> &pos_) const {
    std::vector<size_t> idx;
    for(const auto &p_ : pos_)
        for(const auto &[i, p] : iter::enumerate(pos))
            if(p == p_) { idx.emplace_back(i); }
    return idx;
}

std::vector<size_t> qm::Gate::idx_dn(const std::vector<size_t> &pos_) const {
    std::vector<size_t> idx;
    for(const auto &p_ : pos_)
        for(const auto &[i, p] : iter::enumerate(pos))
            if(p == p_) { idx.emplace_back(i + pos.size()); }
    return idx;
}

qm::Gate qm::Gate::insert(const Gate &other) const { return qm::insert(*this, other); }
qm::Gate qm::Gate::connect_above(const Gate &other) const { return qm::connect(other, *this); }
qm::Gate qm::Gate::connect_under(const Gate &other) const { return qm::connect(*this, other); }

template<auto N>
qm::Gate qm::Gate::trace(const std::array<Eigen::IndexPair<Eigen::Index>, N> &idxpair) const {
    return qm::trace(*this, idxpair);
}
template qm::Gate qm::Gate::trace(const std::array<Eigen::IndexPair<Eigen::Index>, 1> &idxpairs) const;
template qm::Gate qm::Gate::trace(const std::array<Eigen::IndexPair<Eigen::Index>, 2> &idxpairs) const;

qm::Gate qm::Gate::trace_idx(const std::vector<long> &idx_) const { return qm::trace_idx(*this, idx_); }
qm::Gate qm::Gate::trace_pos(const std::vector<size_t> &pos_) const { return qm::trace_pos(*this, pos_); }
qm::Gate qm::Gate::trace_pos(size_t pos_) const { return qm::trace_pos(*this, pos_); }
qm::cplx qm::Gate::trace() const { return qm::trace(*this); }

std::vector<std::vector<size_t>> qm::get_gate_sequence(const std::vector<qm::Gate> &layer) {
    std::vector<std::vector<size_t>> gate_sequence;
    size_t                           gate_size = layer.front().pos.size();
    size_t                           pos_max   = layer.back().pos.back();
    for(size_t offset = 0; offset < gate_size; offset++) {
        if(offset + gate_size > pos_max + 1) break;
        auto off_idx = num::range<size_t>(offset, pos_max - gate_size + 2, gate_size);
        if(num::mod<size_t>(offset, 2) == 1) std::reverse(off_idx.begin(), off_idx.end()); // If odd, reverse the sequence
        gate_sequence.emplace_back(off_idx);
    }
    return gate_sequence;
}

template<iter::order o>
std::vector<std::vector<size_t>> qm::get_lightcone(const std::vector<std::vector<qm::Gate>> &layers, size_t pos) {
    std::vector<std::vector<size_t>> cone;
    cone.emplace_back(std::vector<size_t>{pos});
    for(const auto &[idx_layer, layer] : iter::enumerate<o>(layers)) {
        if(layer.empty()) continue;
        auto gate_sequence = get_gate_sequence(layer);
        for(const auto &[idx_sublayer, seq] : iter::enumerate<o>(gate_sequence)) {
            std::set<size_t> match;
            // The cone width is strictly increasing.
            // So we immediately add the sites from the previous layer.
            match.insert(cone.back().begin(), cone.back().end());
            for(const auto &[idx_seq, pos_gate] : iter::enumerate<o>(seq)) {
                // Now we try to find some more matching sites in the new layer.
                auto               &u = layer[pos_gate];
                std::vector<size_t> pos_isect;
                std::set_intersection(cone.back().begin(), cone.back().end(), u.pos.begin(), u.pos.end(), back_inserter(pos_isect));
                if(not pos_isect.empty()) match.insert(u.pos.begin(), u.pos.end()); // Add positions from u if the gate connects to the current cone
            }
            if(match.empty())
                cone.emplace_back(cone.back());
            else
                cone.emplace_back(std::vector<size_t>(match.begin(), match.end()));
        }
    }
    if constexpr(o == iter::order::rev) std::reverse(cone.begin(), cone.end());
    return cone;
}

template std::vector<std::vector<size_t>> qm::get_lightcone<iter::order::def>(const std::vector<std::vector<qm::Gate>> &layers, size_t pos);
template std::vector<std::vector<size_t>> qm::get_lightcone<iter::order::rev>(const std::vector<std::vector<qm::Gate>> &layers, size_t pos);

std::vector<std::vector<size_t>> qm::get_lightcone_intersection(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t pos_tau, size_t pos_sig) {
    auto tau_cone = qm::get_lightcone<iter::order::def>(unitary_layers, pos_tau);
    auto sig_cone = qm::get_lightcone<iter::order::rev>(unitary_layers, pos_sig); // This cone is upside down!
    if(tau_cone.size() != sig_cone.size()) throw std::runtime_error("tau and sig cones should have equal size!");

    std::vector<std::vector<size_t>> int_cone;

    // Find the intersection between tau and sig cones
    for(size_t idx_sublayer = 0; idx_sublayer < tau_cone.size(); idx_sublayer++) {
        auto               &tau_sublayer = tau_cone[idx_sublayer];
        auto               &sig_sublayer = sig_cone[idx_sublayer];
        std::vector<size_t> pos_isect;
        std::set_intersection(tau_sublayer.begin(), tau_sublayer.end(), sig_sublayer.begin(), sig_sublayer.end(), back_inserter(pos_isect));
        int_cone.emplace_back(pos_isect);
    }
    bool deb = tools::log->level() == spdlog::level::trace;
    if(deb) {
        fmt::print("Lightcones\n");
        auto tau_pic = get_lightcone_picture(unitary_layers, tau_cone, "tau");
        auto sig_pic = get_lightcone_picture(unitary_layers, sig_cone, "sig");
        auto int_pic = get_lightcone_picture(unitary_layers, int_cone, "int");
        for(const auto &c : iter::reverse(tau_pic)) fmt::print("{}\n", c);
        for(const auto &c : iter::reverse(sig_pic)) fmt::print("{}\n", c);
        for(const auto &c : iter::reverse(int_pic)) fmt::print("{}\n", c);
    }

    return int_cone;
}

std::vector<std::string> qm::get_lightcone_picture(const std::vector<std::vector<qm::Gate>> &layers, const std::vector<std::vector<size_t>> &cone,
                                                   std::string_view tag, size_t pw, std::string_view sep) {
    std::vector<std::string> pic;
    if(not layers.empty() and not cone.empty()) {
        size_t sw = sep.size();
        size_t tw = static_cast<size_t>(tag.size()) + 5;                                       // Tag width (brackets, number and colon)
        size_t mw = static_cast<size_t>(layers.front().back().pos.back() + 1) * (pw + sw) + 1; // max cone width
        pic       = std::vector<std::string>(cone.size(), fmt::format("{0:^{1}}", " ", tw + mw));
        for(const auto &[i, c] : iter::enumerate(cone)) {
            pic[i].replace(0, tw, fmt::format("{}[{:^2}]:", tag, i));
            for(const auto &[j, p] : iter::enumerate(c)) pic[i].replace(tw + p * (pw + sw), pw, fmt::format("{0:>{1}}{2}", p, pw, sep));
        }
    }
    return pic;
}

qm::Gate qm::insert(const qm::Gate &middle_gate, const qm::Gate &updown_gate) {
    std::vector<size_t> pos_isect; // locations that intersect: both on middle and updown gates
    std::vector<size_t> pos_nsect; // locations that do not intersect: not in both on middle and updown gates
    std::set_intersection(middle_gate.pos.begin(), middle_gate.pos.end(), updown_gate.pos.begin(), updown_gate.pos.end(), back_inserter(pos_isect));
    std::set_symmetric_difference(middle_gate.pos.begin(), middle_gate.pos.end(), updown_gate.pos.begin(), updown_gate.pos.end(), back_inserter(pos_nsect));
    if(pos_isect.empty()) return middle_gate;
    bool inc = std::includes(middle_gate.pos.begin(), middle_gate.pos.end(), updown_gate.pos.begin(), updown_gate.pos.end());

    auto shp_udn4 = updown_gate.shape<4>();
    auto shp_udn2 = updown_gate.shape<2>();
    if(not pos_isect.empty() and pos_nsect.empty() and inc) {
        // In this case we stack gates vertically that are equally wide.
        // Should be the simplest case
        Eigen::Tensor<cplx, 2> op =
            updown_gate.op.contract(middle_gate.op, tenx::idx({1}, {0})).contract(updown_gate.op.conjugate().shuffle(tenx::array2{1, 0}), tenx::idx({1}, {0}));
        return qm::Gate{op, middle_gate.pos, middle_gate.dim};
    }
    if(pos_isect.size() == 1 and pos_nsect.size() == 1 and middle_gate.pos.size() == 1 and updown_gate.pos.size() == 2) {
        // One common location, one uncommon. Then this connects a 1-site gate with 2-site gates up and down
        tenx::idxlistpair<1> idx1;
        tenx::idxlistpair<2> idx2;
        // Decide if this is connects on the left or right leg
        if(middle_gate.pos.front() == updown_gate.pos.front()) {
            /*  Left insert
             *            0    1              0    1               0    1
             *            |    |              |    |               |    |
             *           [  up  ]            [  up  ]             [  up  ]
             *            |    |              |    |               |    |
             *           (2)   3              |   <2>              |    |
             *           (0)                  |                    |    |
             *            |                   |                    |    |
             *          [mid]       ===>    [mid]        ===>    [mid]  |
             *            |                   |                    |    |
             *            1                  (3)                   |    |
             *            0    1             (0)  <1>              |    |
             *            |    |              |    |               |    |
             *           [  dn  ]            [  dn  ]             [  dn  ]
             *            |    |              |    |               |    |
             *            2    3              2    3               2    3
             *
             */

            idx1 = tenx::idx({2}, {0});
            idx2 = tenx::idx({3, 2}, {2, 3});
        } else {
            /*  Right insert
             *            0    1              0    1               0    1
             *            |    |              |    |               |    |
             *           [  up  ]            [  up  ]             [  up  ]
             *            |    |              |    |               |    |
             *            2   (3)            <2>   |               |    |
             *                (0)                  |               |    |
             *                 |                   |               |    |
             *               [mid]   ===>        [mid]     ===>    |  [mid]
             *                 |                   |               |    |
             *                 1                  (3)              |    |
             *            0    1             <0>  (1)              |    |
             *            |    |              |    |               |    |
             *           [  dn  ]            [  dn  ]             [  dn  ]
             *            |    |              |    |               |    |
             *            2    3              2    3               2    3
             *
             */
            idx1 = tenx::idx({3}, {0});
            idx2 = tenx::idx({2, 3}, {2, 3});
        }
        tools::log->trace("Inserting gate pos {} between gates pos {} | pos_isect {} | pos_nsect {} | inc {} | contract_b", middle_gate.pos, updown_gate.pos,
                          pos_isect, pos_nsect, inc);
        auto op = contract_b(middle_gate.op, updown_gate.op, shp_udn2, shp_udn4, idx1, idx2);
        return qm::Gate{op, updown_gate.pos, updown_gate.dim};
    }
    if(pos_isect.size() == 1 and pos_nsect.size() >= 2 and updown_gate.pos.size() == 2) {
        // One common location, so it must be at the edge.
        std::array<long, 4>         shp_mid4{};
        std::array<long, 2>         dim2{};
        std::array<Eigen::Index, 6> shf6{};
        tenx::idxlistpair<1>        idx1;
        tenx::idxlistpair<2>        idx2;
        std::vector<size_t>         pos;
        std::vector<long>           dim;
        size_t                      merged = pos_nsect.size() - pos_isect.size();
        // Decide if this is connects on the right or left leg
        if(middle_gate.pos.front() == updown_gate.pos.back()) {
            /*  Right insert (Free legs in mid are merged)
             *
             *    0    1                      0    1                        0    1
             *    |    |                      |    |                        |    |
             *   [  up  ]                    [  up  ]                      [  up  ]
             *    |    |                      |    |                        |    |
             *   2    (3)                     |    |                        |    |
             *        (0)    1               <2>   |     3                  |    |     2                                       0   1   2              0
             *         |     |                     |     |                  |    |     |                                       |   |   |              |
             *        [  mid  ]       ===>        [  mid  ]     ===>        |   [  mid  ]       shuffle({0,1,2,4,5,3})  =    [   gate    ]  =   [   gate    ]
             *         |     |                     |     |                  |    |     |                                       |   |   |              |
             *         2     3                    (4)    5                  |    |     3                                       3   4   5              1
             *    0    1                     <0>  (1)                       |    |
             *    |    |                      |    |                        |    |
             *   [  dn  ]                    [  dn  ]                      [  dn  ]
             *    |    |                      |    |                        |    |
             *    2    3                      2    3                        4    5
             *
             */

            if(pos_nsect.size() == 1) shp_mid4 = group(middle_gate.shape<2>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 2) shp_mid4 = group(middle_gate.shape<4>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 3) shp_mid4 = group(middle_gate.shape<6>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 4) shp_mid4 = group(middle_gate.shape<8>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 5) shp_mid4 = group(middle_gate.shape<10>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 6) shp_mid4 = group(middle_gate.shape<12>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 7) shp_mid4 = group(middle_gate.shape<14>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 8) shp_mid4 = group(middle_gate.shape<16>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 9) shp_mid4 = group(middle_gate.shape<18>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 10) shp_mid4 = group(middle_gate.shape<20>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 11) shp_mid4 = group(middle_gate.shape<22>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 12) shp_mid4 = group(middle_gate.shape<24>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 13) shp_mid4 = group(middle_gate.shape<26>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 14) shp_mid4 = group(middle_gate.shape<28>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 15) shp_mid4 = group(middle_gate.shape<30>(), std::array<size_t, 4>{1, merged, 1, merged});
            if(pos_nsect.size() == 16) shp_mid4 = group(middle_gate.shape<32>(), std::array<size_t, 4>{1, merged, 1, merged});
            idx1 = tenx::idx({3}, {0});
            idx2 = tenx::idx({2, 4}, {2, 3});
            shf6 = std::array<Eigen::Index, 6>{0, 1, 2, 4, 5, 3};
            pos  = concat(updown_gate.pos, subset(middle_gate.pos, 1, merged));
            dim  = concat(updown_gate.dim, subset(middle_gate.dim, 1, merged));
            dim2 = repeat(std::array<long, 1>{shp_mid4[1] * shp_udn2[0]});
        } else {
            /*  Left insert (Free legs in mid are merged)
             *
             *           0    1                      0    1                        0    1
             *           |    |                      |    |                        |    |
             *          [  up  ]                    [  up  ]                      [  up  ]
             *           |    |                      |    |                        |    |
             *          (2)   3                      |   <2>                       |    |
             *     0    (1)                    3     |                       2     |    |                                  0   1   2              0
             *     |     |                     |     |                       |     |    |                                  |   |   |              |
             *    [  mid  ]       ===>        [  mid  ]           ===>      [  mid  ]   |   shuffle({2,0,1,3,4,5})  =    [   gate    ]  =   [   gate    ]
             *     |     |                     |     |                       |     |    |                                  |   |   |              |
             *     2     3                     4    (5)                      3     |    |                                  3   4   5              1
             *           0    1                     (0)   <1>                      |    |
             *           |    |                      |    |                        |    |
             *          [  dn  ]                    [  dn  ]                      [  dn  ]
             *           |    |                      |    |                        |    |
             *           2    3                      2    3                        4    5
             *
             */
            if(pos_nsect.size() == 1) shp_mid4 = group(middle_gate.shape<2>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 2) shp_mid4 = group(middle_gate.shape<4>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 3) shp_mid4 = group(middle_gate.shape<6>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 4) shp_mid4 = group(middle_gate.shape<8>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 5) shp_mid4 = group(middle_gate.shape<10>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 6) shp_mid4 = group(middle_gate.shape<12>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 7) shp_mid4 = group(middle_gate.shape<14>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 8) shp_mid4 = group(middle_gate.shape<16>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 9) shp_mid4 = group(middle_gate.shape<18>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 10) shp_mid4 = group(middle_gate.shape<20>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 11) shp_mid4 = group(middle_gate.shape<22>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 12) shp_mid4 = group(middle_gate.shape<24>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 13) shp_mid4 = group(middle_gate.shape<26>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 14) shp_mid4 = group(middle_gate.shape<28>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 15) shp_mid4 = group(middle_gate.shape<30>(), std::array<size_t, 4>{merged, 1, merged, 1});
            if(pos_nsect.size() == 16) shp_mid4 = group(middle_gate.shape<32>(), std::array<size_t, 4>{merged, 1, merged, 1});
            idx1 = tenx::idx({2}, {1});
            idx2 = tenx::idx({5, 2}, {2, 3});

            shf6 = std::array<Eigen::Index, 6>{2, 0, 1, 3, 4, 5};
            pos  = concat(subset(middle_gate.pos, 0, merged), updown_gate.pos);
            dim  = concat(subset(middle_gate.dim, 0, merged), updown_gate.dim);
            dim2 = repeat(std::array<long, 1>{shp_mid4[0] * shp_udn2[0]});
        }
        tools::log->trace("Inserting gate pos {} between gates pos {} | pos_isect {} | pos_nsect {} | inc {} | contract_a", middle_gate.pos, updown_gate.pos,
                          pos_isect, pos_nsect, inc);
        auto op = contract_a(middle_gate.op, updown_gate.op, shp_mid4, shp_udn4, shf6, idx1, idx2, dim2);
        return qm::Gate{op, pos, dim};
    }
    if(pos_isect.size() == 2 and pos_nsect.size() >= 1 and updown_gate.pos.size() == 2) {
        // All the updown locations get contracted
        size_t                      usize = updown_gate.pos.size();
        size_t                      msize = middle_gate.pos.size();
        std::array<long, 4>         shp_mid4{};
        std::array<long, 6>         shp_mid6{};
        std::array<long, 2>         dim2{};
        std::array<Eigen::Index, 6> shf6{};
        std::array<Eigen::Index, 4> shf4{};
        tenx::idxlistpair<1>        idx_up, idx_dn;
        Eigen::Tensor<cplx, 2>      op;
        auto                        offset =
            static_cast<size_t>(std::distance(middle_gate.pos.begin(), find(middle_gate.pos.begin(), middle_gate.pos.end(), updown_gate.pos.front())));
        size_t offmin = 0;
        size_t offmax = middle_gate.pos.size() - updown_gate.pos.size();
        size_t merged = offmax;
        // Decide if this is connects on the left, right or somewhere in the center
        if(offset == offmin or offset == offmax) {
            if(offset == 0) {
                /*  Insert at offset 0
                 *            0                   0                 0
                 *            |                   |                 |
                 *           [up]                [up]              [up]
                 *            |                   |                 |
                 *           (1)                  |                 |
                 *           (0)  1               |   1             |   1
                 *            |   |               |   |             |   |
                 *          [  mid  ]   ===>    [  mid  ]         [  mid  ]     ===> shuffle({0,1,3,2})
                 *           |   |               |   |             |   |
                 *           2   3              (2)   3            |   2
                 *           0                  (0)                |
                 *           |                   |                 |
                 *          [dn]                [dn]              [dn]
                 *           |                   |                 |
                 *           1                   1                 3
                 *
                 */
                if(msize == 2) shp_mid4 = group(middle_gate.shape<4>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 3) shp_mid4 = group(middle_gate.shape<6>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 4) shp_mid4 = group(middle_gate.shape<8>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 5) shp_mid4 = group(middle_gate.shape<10>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 6) shp_mid4 = group(middle_gate.shape<12>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 7) shp_mid4 = group(middle_gate.shape<14>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 8) shp_mid4 = group(middle_gate.shape<16>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 9) shp_mid4 = group(middle_gate.shape<18>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 10) shp_mid4 = group(middle_gate.shape<20>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 11) shp_mid4 = group(middle_gate.shape<22>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 12) shp_mid4 = group(middle_gate.shape<24>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 13) shp_mid4 = group(middle_gate.shape<26>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 14) shp_mid4 = group(middle_gate.shape<28>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 15) shp_mid4 = group(middle_gate.shape<30>(), std::array<size_t, 4>{usize, merged, usize, merged});
                if(msize == 16) shp_mid4 = group(middle_gate.shape<32>(), std::array<size_t, 4>{usize, merged, usize, merged});
                dim2   = middle_gate.shape<2>();
                idx_up = tenx::idx({1}, {0});
                idx_dn = tenx::idx({2}, {1}); // Avoid a shuffle by contracting the other leg. Remember that dn = ud^dagger.

                shf4 = tenx::array4{0, 1, 3, 2};
                tools::log->trace("Inserting gate pos {} between gates pos {} | pos_isect {} | pos_nsect {} | inc {} | contract_d", middle_gate.pos,
                                  updown_gate.pos, pos_isect, pos_nsect, inc);
                op = contract_d(middle_gate.op, updown_gate.op, shp_mid4, idx_up, idx_dn, shf4, dim2);
            } else if(offset == offmax) {
                /*  Insert at offmax
                 *            0                   0                 0
                 *            |                   |                 |
                 *           [up]                [up]              [up]
                 *            |                   |                 |
                 *           (1)                  |                 |
                 *        0  (1)              1   |             1   |
                 *        |   |               |   |             |   |
                 *      [  mid  ]   ===>    [  mid  ]         [  mid  ]     ===> shuffle({1,0,2,3})
                 *        |   |               |   |             |   |
                 *        2   3               2  (3)            2   |
                 *            0                  (0)                |
                 *            |                   |                 |
                 *           [dn]                [dn]              [dn]
                 *            |                   |                 |
                 *            1                   1                 3
                 *
                 */

                if(msize == 2) shp_mid4 = group(middle_gate.shape<4>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 3) shp_mid4 = group(middle_gate.shape<6>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 4) shp_mid4 = group(middle_gate.shape<8>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 5) shp_mid4 = group(middle_gate.shape<10>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 6) shp_mid4 = group(middle_gate.shape<12>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 7) shp_mid4 = group(middle_gate.shape<14>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 8) shp_mid4 = group(middle_gate.shape<16>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 9) shp_mid4 = group(middle_gate.shape<18>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 10) shp_mid4 = group(middle_gate.shape<20>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 11) shp_mid4 = group(middle_gate.shape<22>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 12) shp_mid4 = group(middle_gate.shape<24>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 13) shp_mid4 = group(middle_gate.shape<26>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 14) shp_mid4 = group(middle_gate.shape<28>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 15) shp_mid4 = group(middle_gate.shape<30>(), std::array<size_t, 4>{merged, usize, merged, usize});
                if(msize == 16) shp_mid4 = group(middle_gate.shape<32>(), std::array<size_t, 4>{merged, usize, merged, usize});
                dim2   = middle_gate.shape<2>();
                idx_up = tenx::idx({1}, {1});
                idx_dn = tenx::idx({3}, {1}); // Avoid a shuffle by contracting the other leg. Remember that dn = ud^dagger.
                shf4   = tenx::array4{1, 0, 2, 3};
            }
            tools::log->trace("Inserting gate pos {} between gates pos {} | pos_isect {} | pos_nsect {} | inc {} | contract_d", middle_gate.pos,
                              updown_gate.pos, pos_isect, pos_nsect, inc);
            op = contract_d(middle_gate.op, updown_gate.op, shp_mid4, idx_up, idx_dn, shf4, dim2);
        } else {
            /*  Insert at offmax
             *           0                     0                      0
             *           |                     |                      |
             *          [up]                  [up]                   [up]
             *           |                     |                      |
             *          (1)                    |                      |
             *       0  (1)  2             1   |   2              1   |   2
             *       |   |   |             |   |   |              |   |   |
             *      [   mid   ]    ===>   [   mid   ]    ===>    [   mid   ]       ===> shuffle({1,0,2,3,5,4})
             *       |   |   |             |   |   |              |   |   |
             *       3   4   5             3  (4)  5              3   |   4
             *           0                    (0)                     |
             *           |                     |                      |
             *          [dn]                  [dn]                   [dn]
             *           |                     |                      |
             *           1                     1                      5
             *
             */

            if(msize == 4) shp_mid6 = group(middle_gate.shape<8>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 5) shp_mid6 = group(middle_gate.shape<10>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 6) shp_mid6 = group(middle_gate.shape<12>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 7) shp_mid6 = group(middle_gate.shape<14>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 8) shp_mid6 = group(middle_gate.shape<16>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 9) shp_mid6 = group(middle_gate.shape<18>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 10) shp_mid6 = group(middle_gate.shape<20>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 11) shp_mid6 = group(middle_gate.shape<22>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 12) shp_mid6 = group(middle_gate.shape<24>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 13) shp_mid6 = group(middle_gate.shape<26>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 14) shp_mid6 = group(middle_gate.shape<28>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 15) shp_mid6 = group(middle_gate.shape<30>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            if(msize == 16) shp_mid6 = group(middle_gate.shape<32>(), repeat(std::array<size_t, 3>{offset, usize, msize - usize - offset}));
            dim2   = middle_gate.shape<2>();
            idx_up = tenx::idx({1}, {1});
            idx_dn = tenx::idx({4}, {1}); // Avoid a shuffle by contracting the other leg. Remember that dn = ud^dagger.
            shf6   = tenx::array6{1, 0, 2, 3, 5, 4};
            tools::log->trace("Inserting gate pos {} between gates pos {} | pos_isect {} | pos_nsect {} | inc {} | contract_c", middle_gate.pos,
                              updown_gate.pos, pos_isect, pos_nsect, inc);
            op = contract_c(middle_gate.op, updown_gate.op, shp_mid6, idx_up, idx_dn, shf6, dim2);
        }
        if(op.size() == 0) throw std::logic_error("No op computed!");
        return qm::Gate{op, middle_gate.pos, middle_gate.dim};
    }

    throw std::runtime_error(fmt::format("Insert case not implemented: middle pos {} | updown pos {} | pos_isect {} | pos_nsect {}", middle_gate.pos,
                                         updown_gate.pos, pos_isect, pos_nsect));
}

qm::Gate qm::connect(const qm::Gate &dn_gate, const qm::Gate &up_gate) {
    std::vector<size_t> pos_isect; // locations that intersect: both on middle and updown gates
    std::vector<size_t> pos_nsect; // locations that do not intersect: not in both on middle and updown gates
    std::set_intersection(dn_gate.pos.begin(), dn_gate.pos.end(), up_gate.pos.begin(), up_gate.pos.end(), back_inserter(pos_isect));
    std::set_symmetric_difference(dn_gate.pos.begin(), dn_gate.pos.end(), up_gate.pos.begin(), up_gate.pos.end(), back_inserter(pos_nsect));
    bool inc = std::includes(dn_gate.pos.begin(), dn_gate.pos.end(), up_gate.pos.begin(), up_gate.pos.end());
    if(pos_isect.empty()) return dn_gate;

    tools::log->trace("Connecting dn {} | up {}", dn_gate.pos, up_gate.pos);

    if(not pos_isect.empty() and pos_nsect.empty() and inc) {
        // This case is a vertical stack. Should be the simplest case
        Eigen::Tensor<cplx, 2> op = up_gate.op.contract(dn_gate.op, tenx::idx({1}, {0}));
        return qm::Gate{op, dn_gate.pos, dn_gate.dim};
    }
    if(not pos_isect.empty() and not pos_nsect.empty() and inc) {
        std::array<long, 6> dn_shp6{};
        std::array<long, 4> dn_shp4{};
        std::vector<size_t> pos;
        std::vector<long>   dim;
        auto                dn_size  = dn_gate.pos.size();
        auto                up_size  = up_gate.pos.size();
        size_t              dn_merge = dn_size - pos_isect.size();
        size_t              up_merge = up_size;
        size_t              offset   = up_gate.pos.front() - dn_gate.pos.front();
        // Decide if this is connects on the right or right leg
        if(dn_gate.pos.back() == pos_isect.back()) {
            /*  Right connection
             *
             *            |    |
             *            |  [up]
             *            |   |
             *          [  dn  ]
             *           |    |
             *
             */
            if(dn_size == 2) dn_shp4 = group(dn_gate.shape<4>(), repeat(std::array<size_t, 2>{dn_merge, up_merge}));
            if(dn_size == 3) dn_shp4 = group(dn_gate.shape<6>(), repeat(std::array<size_t, 2>{dn_merge, up_merge}));
            if(dn_size == 4) dn_shp4 = group(dn_gate.shape<8>(), repeat(std::array<size_t, 2>{dn_merge, up_merge}));
            if(dn_size == 5) dn_shp4 = group(dn_gate.shape<10>(), repeat(std::array<size_t, 2>{dn_merge, up_merge}));
            if(dn_size == 6) dn_shp4 = group(dn_gate.shape<12>(), repeat(std::array<size_t, 2>{dn_merge, up_merge}));
            if(dn_size == 7) dn_shp4 = group(dn_gate.shape<14>(), repeat(std::array<size_t, 2>{dn_merge, up_merge}));
            if(dn_size == 8) dn_shp4 = group(dn_gate.shape<16>(), repeat(std::array<size_t, 2>{dn_merge, up_merge}));
            auto                   dim2 = dn_gate.shape<2>();
            Eigen::Tensor<cplx, 2> op   = up_gate.op.contract(dn_gate.op.reshape(dn_shp4), tenx::idx({1}, {1})).shuffle(tenx::array4{1, 0, 2, 3}).reshape(dim2);
            return qm::Gate{op, dn_gate.pos, dn_gate.dim};
        } else if(dn_gate.pos.front() == pos_isect.front()) {
            /*  Left connection
             *
             *            |   |
             *          [up]  |
             *           |    |
             *          [  dn  ]
             *           |    |
             */
            if(dn_size == 2) dn_shp4 = group(dn_gate.shape<4>(), repeat(std::array<size_t, 2>{up_merge, dn_merge}));
            if(dn_size == 3) dn_shp4 = group(dn_gate.shape<6>(), repeat(std::array<size_t, 2>{up_merge, dn_merge}));
            if(dn_size == 4) dn_shp4 = group(dn_gate.shape<8>(), repeat(std::array<size_t, 2>{up_merge, dn_merge}));
            if(dn_size == 5) dn_shp4 = group(dn_gate.shape<10>(), repeat(std::array<size_t, 2>{up_merge, dn_merge}));
            if(dn_size == 6) dn_shp4 = group(dn_gate.shape<12>(), repeat(std::array<size_t, 2>{up_merge, dn_merge}));
            if(dn_size == 7) dn_shp4 = group(dn_gate.shape<14>(), repeat(std::array<size_t, 2>{up_merge, dn_merge}));
            if(dn_size == 8) dn_shp4 = group(dn_gate.shape<16>(), repeat(std::array<size_t, 2>{up_merge, dn_merge}));
            auto                   dim2 = dn_gate.shape<2>();
            Eigen::Tensor<cplx, 2> op   = up_gate.op.contract(dn_gate.op.reshape(dn_shp4), tenx::idx({1}, {0})).reshape(dim2);
            return qm::Gate{op, dn_gate.pos, dn_gate.dim};
        } else {
            /*  Center connection
             *
             *           |    |   |
             *           |  [up]  |
             *           |   |    |
             *          [    dn    ]
             *           |    |   |
             */
            if(dn_size == 3) dn_shp6 = group(dn_gate.shape<6>(), repeat(std::array<size_t, 3>{offset, up_size, dn_size - up_size - offset}));
            if(dn_size == 4) dn_shp6 = group(dn_gate.shape<8>(), repeat(std::array<size_t, 3>{offset, up_size, dn_size - up_size - offset}));
            if(dn_size == 5) dn_shp6 = group(dn_gate.shape<10>(), repeat(std::array<size_t, 3>{offset, up_size, dn_size - up_size - offset}));
            if(dn_size == 6) dn_shp6 = group(dn_gate.shape<12>(), repeat(std::array<size_t, 3>{offset, up_size, dn_size - up_size - offset}));
            if(dn_size == 7) dn_shp6 = group(dn_gate.shape<14>(), repeat(std::array<size_t, 3>{offset, up_size, dn_size - up_size - offset}));
            if(dn_size == 8) dn_shp6 = group(dn_gate.shape<16>(), repeat(std::array<size_t, 3>{offset, up_size, dn_size - up_size - offset}));
            auto                   dim2 = dn_gate.shape<2>();
            Eigen::Tensor<cplx, 2> op =
                up_gate.op.contract(dn_gate.op.reshape(dn_shp6), tenx::idx({1}, {1})).shuffle(tenx::array6{1, 0, 2, 3, 4, 5}).reshape(dim2);
            return qm::Gate{op, dn_gate.pos, dn_gate.dim};
        }
    }
    if(pos_isect.size() == 1 and pos_nsect.size() >= 2) {
        // One common location, and one uncommon. Then this connects two 2-site gates
        std::array<long, 4> dn_shp4{};
        std::array<long, 4> up_shp4{};
        std::vector<size_t> pos;
        std::vector<long>   dim;
        auto                dn_size  = dn_gate.pos.size();
        auto                up_size  = up_gate.pos.size();
        size_t              dn_merge = dn_size - pos_isect.size();
        size_t              up_merge = up_size - pos_isect.size();
        // Decide if this is connects on the right or right leg
        if(dn_gate.pos.front() == up_gate.pos.back()) {
            /*  Right connection
             *
             *      |    |     |           0   1   2              0
             *     [  up  ]    |           |   |   |              |
             *      |    |     |   =     [   gate    ]  =   [   gate    ]
             *      |   [  dn  ]           |   |   |              |
             *      |    |     |           3   4   5              1
             *
             */
            if(dn_size == 2) dn_shp4 = group(dn_gate.shape<4>(), repeat(std::array<size_t, 2>{1, dn_merge}));
            if(dn_size == 3) dn_shp4 = group(dn_gate.shape<6>(), repeat(std::array<size_t, 2>{1, dn_merge}));
            if(dn_size == 4) dn_shp4 = group(dn_gate.shape<8>(), repeat(std::array<size_t, 2>{1, dn_merge}));
            if(dn_size == 5) dn_shp4 = group(dn_gate.shape<10>(), repeat(std::array<size_t, 2>{1, dn_merge}));
            if(dn_size == 6) dn_shp4 = group(dn_gate.shape<12>(), repeat(std::array<size_t, 2>{1, dn_merge}));
            if(dn_size == 7) dn_shp4 = group(dn_gate.shape<14>(), repeat(std::array<size_t, 2>{1, dn_merge}));
            if(dn_size == 8) dn_shp4 = group(dn_gate.shape<16>(), repeat(std::array<size_t, 2>{1, dn_merge}));
            if(up_size == 2) up_shp4 = group(dn_gate.shape<4>(), repeat(std::array<size_t, 2>{up_merge, 1}));
            if(up_size == 3) up_shp4 = group(dn_gate.shape<6>(), repeat(std::array<size_t, 2>{up_merge, 1}));
            if(up_size == 4) up_shp4 = group(dn_gate.shape<8>(), repeat(std::array<size_t, 2>{up_merge, 1}));
            if(up_size == 5) up_shp4 = group(dn_gate.shape<10>(), repeat(std::array<size_t, 2>{up_merge, 1}));
            if(up_size == 6) up_shp4 = group(dn_gate.shape<12>(), repeat(std::array<size_t, 2>{up_merge, 1}));
            if(up_size == 7) up_shp4 = group(dn_gate.shape<14>(), repeat(std::array<size_t, 2>{up_merge, 1}));
            if(up_size == 8) up_shp4 = group(dn_gate.shape<16>(), repeat(std::array<size_t, 2>{up_merge, 1}));
            auto                   dim2 = repeat(std::array<long, 1>{up_shp4[0] * up_shp4[1] * dn_shp4[1]});
            Eigen::Tensor<cplx, 2> op =
                dn_gate.op.reshape(dn_shp4).contract(up_gate.op.reshape(up_shp4), tenx::idx({0}, {3})).shuffle(tenx::array6{3, 4, 0, 5, 1, 2}).reshape(dim2);
            pos = concat(up_gate.pos, subset(dn_gate.pos, 1, dn_merge));
            dim = concat(up_gate.dim, subset(dn_gate.dim, 1, dn_merge));
            return qm::Gate{op, pos, dim};
        } else {
            /*  Left connection
             *
             *     |     |    |             0   1   2              0
             *     |    [  up  ]            |   |   |              |
             *     |     |    |      =    [   gate    ]  =   [   gate    ]
             *    [  dn  ]    |             |   |   |              |
             *     |     |    |             3   4   5              1
             */
            if(dn_size == 2) dn_shp4 = group(dn_gate.shape<4>(), repeat(std::array<size_t, 2>{dn_merge, 1}));
            if(dn_size == 3) dn_shp4 = group(dn_gate.shape<6>(), repeat(std::array<size_t, 2>{dn_merge, 1}));
            if(dn_size == 4) dn_shp4 = group(dn_gate.shape<8>(), repeat(std::array<size_t, 2>{dn_merge, 1}));
            if(dn_size == 5) dn_shp4 = group(dn_gate.shape<10>(), repeat(std::array<size_t, 2>{dn_merge, 1}));
            if(dn_size == 6) dn_shp4 = group(dn_gate.shape<12>(), repeat(std::array<size_t, 2>{dn_merge, 1}));
            if(dn_size == 7) dn_shp4 = group(dn_gate.shape<14>(), repeat(std::array<size_t, 2>{dn_merge, 1}));
            if(dn_size == 8) dn_shp4 = group(dn_gate.shape<16>(), repeat(std::array<size_t, 2>{dn_merge, 1}));
            if(up_size == 2) up_shp4 = group(dn_gate.shape<4>(), repeat(std::array<size_t, 2>{1, up_merge}));
            if(up_size == 3) up_shp4 = group(dn_gate.shape<6>(), repeat(std::array<size_t, 2>{1, up_merge}));
            if(up_size == 4) up_shp4 = group(dn_gate.shape<8>(), repeat(std::array<size_t, 2>{1, up_merge}));
            if(up_size == 5) up_shp4 = group(dn_gate.shape<10>(), repeat(std::array<size_t, 2>{1, up_merge}));
            if(up_size == 6) up_shp4 = group(dn_gate.shape<12>(), repeat(std::array<size_t, 2>{1, up_merge}));
            if(up_size == 7) up_shp4 = group(dn_gate.shape<14>(), repeat(std::array<size_t, 2>{1, up_merge}));
            if(up_size == 8) up_shp4 = group(dn_gate.shape<16>(), repeat(std::array<size_t, 2>{1, up_merge}));
            auto                   dim2 = repeat(std::array<long, 1>{dn_shp4[0] * up_shp4[0] * up_shp4[1]});
            Eigen::Tensor<cplx, 2> op =
                dn_gate.op.reshape(dn_shp4).contract(up_gate.op.reshape(up_shp4), tenx::idx({1}, {2})).shuffle(tenx::array6{0, 3, 4, 1, 2, 5}).reshape(dim2);
            pos = concat(subset(dn_gate.pos, 0, dn_merge), up_gate.pos);
            dim = concat(subset(dn_gate.dim, 0, dn_merge), up_gate.dim);
            return qm::Gate{op, pos, dim};
        }
    }

    throw std::runtime_error(
        fmt::format("Connect case not implemented: dn pos {} | up pos {} | pos_isect {} | pos_nsect {}", dn_gate.pos, up_gate.pos, pos_isect, pos_nsect));
}

qm::Gate qm::trace_idx(const qm::Gate &gate, const std::vector<long> &idx) {
    if(idx.size() == 2) return qm::trace(gate, tenx::idx({idx[0]}, {idx[1]}));
    // idx is already pair-ordered, so when we call tenx::idx we need to undo that order
    // which explains the weird indexing in the following case
    if(idx.size() == 4) return qm::trace(gate, tenx::idx({idx[0], idx[2]}, {idx[1], idx[3]}));
    throw std::runtime_error(fmt::format("Tracing {} indices is not implemented", idx.size()));
}

qm::Gate qm::trace_pos(const qm::Gate &gate, const std::vector<size_t> &pos) {
    if(pos.size() <= 2) {
        auto              idx = gate.idx(pos);
        std::vector<long> idx_long(idx.begin(), idx.end());
        return qm::trace_idx(gate, idx_long);
    }
    qm::Gate tmp = gate;
    for(auto &p : pos) { tmp = tmp.trace_pos(std::vector<size_t>{p}); }
    return tmp;
}

qm::Gate qm::trace_pos(const qm::Gate &gate, size_t pos) { return qm::trace_pos(gate, gate.idx(std::vector<size_t>{pos})); }

qm::cplx qm::trace(const qm::Gate &gate) {
    qm::Gate t = qm::trace_pos(gate, gate.pos);
    if(not t.pos.empty()) throw std::logic_error(fmt::format("Gate should be empty after tracing all positions. Got pos: {}", t.pos));
    if(t.op.dimension(0) * t.op.dimension(1) != 1)
        throw std::logic_error(fmt::format("Gate should have cplx op after tracing all positions. Got dims: {}", t.op.dimensions()));
    return t.op(0);
}

template<auto N>
qm::Gate qm::trace(const qm::Gate &gate, const std::array<Eigen::IndexPair<Eigen::Index>, N> &idxpairs) {
    // Trace off pairs of indices on a gate
    // This is done one by one recursively until all pairs have been traced
    if constexpr(N == 1) tools::log->trace("Tracing pos {} | dim {} | index pair [{},{}]", gate.pos, gate.dim, idxpairs[0].first, idxpairs[0].second);
    if constexpr(N == 2)
        tools::log->trace("Tracing gate pos {} | dim {} | index pairs [{},{}][{},{}]", gate.pos, gate.dim, idxpairs[0].first, idxpairs[0].second,
                          idxpairs[1].first, idxpairs[1].second);

    // Compute the remaining indices positions and dimensions
    auto idx = gate.idx(); // Twice as long as pos and dim! It has the top and bottom indices
    auto pos = gate.pos;
    auto dim = gate.dim;
    for(const auto &pair : idxpairs) { // Idx pairs has the indices of the legs we remove
        auto it1 = std::find(idx.begin(), idx.end(), std::max(pair.first, pair.second));
        if(it1 != idx.end()) {
            auto idx_rm1 = std::distance(idx.begin(), it1);
            if(idx_rm1 < static_cast<long>(pos.size())) pos.erase(pos.begin() + idx_rm1);
            if(idx_rm1 < static_cast<long>(dim.size())) dim.erase(dim.begin() + idx_rm1);
            idx.erase(it1);
        }
        auto it2 = std::find(idx.begin(), idx.end(), std::min(pair.first, pair.second));
        if(it2 != idx.end()) {
            auto idx_rm2 = std::distance(idx.begin(), it2);
            if(idx_rm2 < static_cast<long>(pos.size())) pos.erase(pos.begin() + idx_rm2);
            if(idx_rm2 < static_cast<long>(dim.size())) dim.erase(dim.begin() + idx_rm2);
            idx.erase(it2);
        }
    }
    if(idx.size() != 2 * (gate.pos.size() - N)) throw std::logic_error("Wrong number of indices removed");

    // Assuming the gate is symmetric, we can compute the rank2 dimensions for the storage of the gate op
    long dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies<>()); // Product of all dimensions of the remaining top legs of the gate
    std::array<long, 2> dim2{dim_prod, dim_prod};

    // Shorthand tensors
    using T2  = Eigen::Tensor<cplx, 2>;
    using T4  = Eigen::Tensor<cplx, 4>;
    using T6  = Eigen::Tensor<cplx, 6>;
    using T8  = Eigen::Tensor<cplx, 8>;
    using T10 = Eigen::Tensor<cplx, 10>;
    using T12 = Eigen::Tensor<cplx, 12>;
    using T14 = Eigen::Tensor<cplx, 14>;
    using T16 = Eigen::Tensor<cplx, 16>;
    using T18 = Eigen::Tensor<cplx, 18>;
    using T20 = Eigen::Tensor<cplx, 20>;
    using T22 = Eigen::Tensor<cplx, 22>;
    using T24 = Eigen::Tensor<cplx, 24>;
    using T26 = Eigen::Tensor<cplx, 26>;
    using T28 = Eigen::Tensor<cplx, 28>;
    using T30 = Eigen::Tensor<cplx, 30>;
    using T32 = Eigen::Tensor<cplx, 32>;
    Eigen::Tensor<cplx, 2> op_traced;

    // Trace
    if(gate.dim.size() == 1) {
        if constexpr(N == 1) { // dim size 1 is special! It can't take 2 index pairs
            op_traced = linalg::tensor::trace(static_cast<T2>(gate.op.reshape(gate.shape<2>())), idxpairs).reshape(dim2);
        } else {
            throw std::runtime_error(fmt::format("Can't trace {} index pairs on gate with pos {}", N, gate.pos));
        }
    }
    /* clang-format off */
    else if(gate.dim.size() == 2 ) op_traced = linalg::tensor::trace(static_cast<T4>(gate.op.reshape(gate.shape<4>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 3 ) op_traced = linalg::tensor::trace(static_cast<T6>(gate.op.reshape(gate.shape<6>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 4 ) op_traced = linalg::tensor::trace(static_cast<T8>(gate.op.reshape(gate.shape<8>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 5 ) op_traced = linalg::tensor::trace(static_cast<T10>(gate.op.reshape(gate.shape<10>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 6 ) op_traced = linalg::tensor::trace(static_cast<T12>(gate.op.reshape(gate.shape<12>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 7 ) op_traced = linalg::tensor::trace(static_cast<T14>(gate.op.reshape(gate.shape<14>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 8 ) op_traced = linalg::tensor::trace(static_cast<T16>(gate.op.reshape(gate.shape<16>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 9 ) op_traced = linalg::tensor::trace(static_cast<T18>(gate.op.reshape(gate.shape<18>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 10) op_traced = linalg::tensor::trace(static_cast<T20>(gate.op.reshape(gate.shape<20>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 11) op_traced = linalg::tensor::trace(static_cast<T22>(gate.op.reshape(gate.shape<22>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 12) op_traced = linalg::tensor::trace(static_cast<T24>(gate.op.reshape(gate.shape<24>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 13) op_traced = linalg::tensor::trace(static_cast<T26>(gate.op.reshape(gate.shape<26>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 14) op_traced = linalg::tensor::trace(static_cast<T28>(gate.op.reshape(gate.shape<28>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 15) op_traced = linalg::tensor::trace(static_cast<T30>(gate.op.reshape(gate.shape<30>())), idxpairs).reshape(dim2);
    else if(gate.dim.size() == 16) op_traced = linalg::tensor::trace(static_cast<T32>(gate.op.reshape(gate.shape<32>())), idxpairs).reshape(dim2);
    /* clang-format on */
    else
        throw std::runtime_error(fmt::format("Trace not implemented: N == {} | dim.size() == {}", N, gate.dim.size()));

    return Gate{op_traced, pos, dim};
}

template qm::Gate qm::trace(const qm::Gate &gate, const std::array<Eigen::IndexPair<Eigen::Index>, 1> &idxpairs);
template qm::Gate qm::trace(const qm::Gate &gate, const std::array<Eigen::IndexPair<Eigen::Index>, 2> &idxpairs);
