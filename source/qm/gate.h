#pragma once
#include "debug/exceptions.h"
#include "math/float.h"
#include "math/tenx.h"
#include "qm.h"
#include <array>
#include <complex>
#include <deque>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

namespace iter {
    enum class order;
}
namespace Eigen {
    template<typename Idx>
    struct IndexPair;
}
enum class GateOp;

namespace qm {

    /* clang-format off */
    class Gate;
    template<typename L>
    [[nodiscard]] std::vector<std::vector<size_t>> get_gate_sequence(const L & layer, bool reverse_odd = true);
    template<iter::order o>
    [[nodiscard]] std::vector<std::vector<size_t>>   get_lightcone(const std::vector<std::vector<qm::Gate>> & layers, size_t pos);
    [[nodiscard]] std::vector<std::vector<size_t>>   get_lightcone_intersection(const std::vector<std::vector<qm::Gate>> & layers, size_t pos_tau, size_t pos_sig);
    [[nodiscard]] std::vector<std::deque<qm::Gate>>  get_lightcone_gate_selection(const std::vector<std::vector<qm::Gate>> & layers,
                                                                                  const std::vector<std::vector<size_t>> & lightcone_intersection, bool reverse_odd = true);
    [[nodiscard]] std::vector<std::deque<qm::Gate>> get_lightcone_gate_selection(const std::vector<std::vector<qm::Gate>> & layers, size_t pos_tau, size_t pos_sig, bool reverse_odd = true);
    [[nodiscard]] std::vector<std::string> get_lightcone_picture(const std::vector<std::vector<qm::Gate>> & layers,
                                                                const std::vector<std::vector<size_t>> & cone, std::string_view tag,
                                                                size_t point_width = 3, std::string_view sep = ",");
    [[nodiscard]] qm::Gate insert(const qm::Gate & middle_gate , const qm::Gate & updown_gate);
    [[nodiscard]] qm::Gate connect(const qm::Gate & dn_gate , const qm::Gate & up_gate);
    [[nodiscard]] qm::Gate trace_idx(const qm::Gate & gate , const std::vector<long> & idx);
    [[nodiscard]] qm::Gate trace_pos(const qm::Gate & gate , const std::vector<size_t> & pos);
    [[nodiscard]] qm::Gate trace_pos(const qm::Gate & gate , size_t pos);
    [[nodiscard]] cplx   trace(const qm::Gate & gate);
    template<auto N>
    [[nodiscard]] qm::Gate trace(const qm::Gate & gate , const std::array<Eigen::IndexPair<Eigen::Index>, N>  & idxpair);


    class Gate {
        public:
        protected:
        template<typename scalar_t, typename alpha_t>
        Eigen::Tensor<scalar_t, 2> exp_internal(const Eigen::Tensor<scalar_t,2> & op_, alpha_t alpha) const;

        mutable std::optional<std::vector<Eigen::Tensor<cplx,2>>> op_split;
        mutable std::optional<Eigen::Tensor<cplx,2>> cnj = std::nullopt;
        mutable std::optional<Eigen::Tensor<cplx,2>> adj = std::nullopt;
        mutable std::optional<Eigen::Tensor<cplx,2>> trn = std::nullopt;
        mutable bool used = false;
        public:
        Eigen::Tensor<cplx,2> op;
        Eigen::Tensor<cplx_t,2> op_t;
        std::vector<size_t> pos;
        std::vector<long> dim;
        Gate() = default;

        template<typename T, typename Device = Eigen::DefaultDevice>
        Gate(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &op_, std::vector<size_t> pos_, std::vector<long> dim_, const Device &device = Device())
            : pos(std::move(pos_)), dim(std::move(dim_)) {
            auto tensor   = tenx::asEval(op_, device);
            auto dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies<>());
            static_assert(tensor.rank() == 2);
            if(dim_prod != tensor.dimension(0) or dim_prod != tensor.dimension(1))
                throw except::logic_error("dim {} not compatible with matrix dimensions {} x {}", dim, tensor.dimension(0), tensor.dimension(1));
            if(pos.size() != dim.size()) throw except::logic_error("pos.size() {} != dim.size() {}", pos, dim);
            // We use a unary expression to cast from std::complex<__float128> to std::complex<double>
            op = op_.unaryExpr([](auto z){return std::complex<real>(static_cast<real>(z.real()), static_cast<real>(z.imag()));});
            op_t = op_.template cast<cplx_t>();
        }
        template<typename T>
        Gate(const Eigen::EigenBase<T> & op_, std::vector<size_t> pos_, std::vector<long> dim_)
            : Gate(tenx::TensorCast(op_), std::move(pos_), std::move(dim_)) {}

        explicit Gate(const Eigen::Tensor<cplx,2> & op_, std::vector<size_t> pos_, std::vector<long> dim_, cplx alpha);
        explicit Gate(const Eigen::Tensor<cplx,2> & op_, std::vector<size_t> pos_, std::vector<long> dim_, cplx_t alpha);
        explicit Gate(const Eigen::Tensor<cplx_t,2> & op_, std::vector<size_t> pos_, std::vector<long> dim_, cplx alpha);
        explicit Gate(const Eigen::Tensor<cplx_t,2> & op_, std::vector<size_t> pos_, std::vector<long> dim_, cplx_t alpha);

        void mark_as_used() const;
        void unmark_as_used() const;
        bool was_used() const;
        [[nodiscard]] Gate exp(cplx alpha) const;
        [[nodiscard]] Gate exp(cplx_t alpha) const;
        [[nodiscard]] bool isUnitary(double prec = 1e-12) const;
        [[nodiscard]] const Eigen::Tensor<cplx,2>& conjugate() const;
        [[nodiscard]] const Eigen::Tensor<cplx,2>& transpose() const;
        [[nodiscard]] const Eigen::Tensor<cplx,2>& adjoint() const;
        [[nodiscard]] const Eigen::Tensor<cplx,2>& unaryOp(GateOp unop) const;
        [[nodiscard]] Gate insert(const Gate & other) const;
        [[nodiscard]] Gate connect_above(const Gate & other) const;
        [[nodiscard]] Gate connect_below(const Gate & other) const;
        template<auto N>
        [[nodiscard]] Gate trace(const std::array<Eigen::IndexPair<Eigen::Index>, N> & idxpair) const;
        [[nodiscard]] Gate trace_idx(const std::vector<long> & idx_) const;
        [[nodiscard]] Gate trace_pos(const std::vector<size_t> & pos_) const;
        [[nodiscard]] Gate trace_pos(size_t pos_) const;
        [[nodiscard]] cplx trace() const;
        template<auto rank>
        [[nodiscard]] std::array<long,rank> shape() const;
        [[nodiscard]] std::vector<size_t> idx() const;
        [[nodiscard]] std::vector<size_t> idx_up() const;
        [[nodiscard]] std::vector<size_t> idx_dn() const;
        [[nodiscard]] std::vector<size_t> idx(const std::vector<size_t> &pos_) const;
        [[nodiscard]] std::vector<size_t> idx_up(const std::vector<size_t> &pos_) const;
        [[nodiscard]] std::vector<size_t> idx_dn(const std::vector<size_t> &pos_) const;
        [[nodiscard]] bool has_pos(const std::vector<size_t> &pos_) const;
        [[nodiscard]] std::vector<size_t> pos_intersection(const std::vector<size_t> &pos_) const;
        [[nodiscard]] std::vector<size_t> pos_difference(const std::vector<size_t> &pos_) const;
        void draw_pos(std::string & layer_str, std::optional<std::string> layer_tag = std::nullopt) const;
    };
    /* clang-format on */
    struct Rwap;
    struct Swap {
        size_t posL, posR;
        Swap() = default;
        Swap(size_t posL, size_t posR);
        bool operator==(const Rwap &rwap) const;
    };
    struct Rwap {
        size_t posL, posR;
        Rwap() = default;
        Rwap(size_t posL, size_t posR);
        bool operator==(const Swap &swap) const;
    };
    class SwapGate : public Gate {
        using Gate::Gate;

        public:
        std::deque<Swap>       swaps;
        std::deque<Rwap>       rwaps; // swaps sequence and reverse swap sequence
        [[nodiscard]] SwapGate exp(cplx alpha) const;
        [[nodiscard]] SwapGate exp(cplx_t alpha) const;
        void                   generate_swap_sequences();
        size_t                 cancel_swaps(std::deque<Rwap> &other_rwaps);
        size_t                 cancel_rwaps(std::deque<Swap> &other_swaps);
    };

    class MpoGate : public Gate {
        protected:
        std::vector<Eigen::Tensor<cplx, 4>> ops;

        public:
        [[nodiscard]] MpoGate                insert(const MpoGate &other) const;
        [[nodiscard]] MpoGate                insert(const Gate &other) const;
        [[nodiscard]] Eigen::Tensor<cplx, 4> get_op();
        [[nodiscard]] Eigen::Tensor<cplx, 4> get_op(const std::vector<size_t> &pos);
        [[nodiscard]] Gate                   split() const; /*!< Split into individual sites with svd */
    };

}