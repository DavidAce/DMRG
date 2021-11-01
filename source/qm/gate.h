#pragma once
#include "qm.h"
#include <array>
#include <complex>
#include <deque>
#include <math/tenx/fwd_decl.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

namespace iter {
    enum class order;
}
namespace Eigen {
    template<typename Idx>
    struct IndexPair;
}

namespace qm {

    /* clang-format off */
    class Gate;
    [[nodiscard]] std::vector<std::vector<size_t>> get_gate_sequence(const std::vector<qm::Gate> & layer);
    template<iter::order o>
    [[nodiscard]] std::vector<std::vector<size_t>> get_lightcone(const std::vector<std::vector<qm::Gate>> & layers, size_t pos);
    [[nodiscard]] std::vector<std::vector<size_t>> get_lightcone_intersection(const std::vector<std::vector<qm::Gate>> & layers, size_t pos_tau, size_t pos_sig);
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
        protected:
        Eigen::Tensor<cplx,2> exp_internal(const Eigen::Tensor<cplx,2> & op_, cplx alpha) const;
        mutable std::optional<Eigen::Tensor<cplx,2>> adj = std::nullopt;
        mutable bool used = false;
        public:
        Eigen::Tensor<cplx,2> op;
        std::vector<size_t> pos;
        std::vector<long> dim;
        Gate() = default;
        Gate(const Eigen::Matrix<cplx,Eigen::Dynamic,Eigen::Dynamic> & op_, std::vector<size_t> pos_, std::vector<long> dim_);
        Gate(const Eigen::Tensor<cplx,2> & op_, std::vector<size_t> pos_, std::vector<long> dim_);
        Gate(const Eigen::Tensor<cplx,2> & op_, std::vector<size_t> pos_, std::vector<long> dim_, cplx alpha);
        void exp_inplace(cplx alpha);
        void mark_as_used() const;
        void unmark_as_used() const;
        bool was_used() const;
        [[nodiscard]] Gate exp(cplx alpha) const;
        [[nodiscard]] bool isUnitary(double prec = 1e-12) const;
        [[nodiscard]] Eigen::Tensor<cplx,2>& adjoint() const;
        [[nodiscard]] Gate insert(const Gate & other) const;
        [[nodiscard]] Gate connect_above(const Gate & other) const;
        [[nodiscard]] Gate connect_under(const Gate & other) const;
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
    };
    /* clang-format on */
    struct Rwap;
    struct Swap {
        size_t posL, posR;
        Swap(size_t posL, size_t posR);
        bool operator==(const Rwap &rwap) const;
    };
    struct Rwap {
        size_t posL, posR;
        Rwap(size_t posL, size_t posR);
        bool operator==(const Swap &swap) const;
    };
    class SwapGate : public Gate {
        using Gate::Gate;

        public:
        std::deque<Swap>       swaps;
        std::deque<Rwap>       rwaps; // swaps sequence and reverse swap sequence
        [[nodiscard]] SwapGate exp(cplx alpha) const;
        void                   generate_swap_sequences();
        long                   cancel_swaps(std::deque<Rwap> &other_rwaps);
        long                   cancel_rwaps(std::deque<Swap> &other_swaps);
    };

}