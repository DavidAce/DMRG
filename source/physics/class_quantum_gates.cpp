//
// Created by david on 2020-12-20.
//

#include "class_quantum_gates.h"
#include "general/nmspc_tensor_extra.h"
#include <Eigen/Core>
#include <general/nmspc_iter.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <utility>

template<typename T>
std::vector<T> subset(const std::vector<T> & vec, size_t idx_start, size_t num){
    if(idx_start + num > vec.size()) throw std::range_error(fmt::format("Vector subset start {} num {} out of range for vector of size {}", idx_start, num, vec.size()));
    auto vec_bgn = vec.begin() + static_cast<long>(idx_start);
    auto vec_end   = vec_bgn + static_cast<long>(num);
    return std::vector<T>(vec_bgn, vec_end);
}


template<auto N, typename T,auto M>
std::array<T,N> subset(const std::array<T,M> &arr, size_t idx_start){
    if(idx_start + N > arr.size()) throw std::range_error(fmt::format("Vector subset start {} num {} out of range for vector of size {}", idx_start,  N, arr.size()));
    auto vec_bgn = arr.begin() + idx_start;
    auto vec_end   = vec_bgn +  N;
    std::array<T,N> res{};
    std::copy(vec_bgn,vec_end, res.begin());
    return res;
}


template<typename T>
T prod(const std::vector<T> & vec, size_t idx_start, size_t num, T seed = 1){
    auto sub = subset(vec, idx_start, num);
    return std::accumulate(sub.begin(),sub.end(), std::multiplies(seed));
}

template<auto N, typename T, auto M>
T prod(const std::array<T,M> & arr, size_t idx_start, T seed = 1){
    auto sub = subset<N>(arr, idx_start);
    T acc = seed;
    for(auto & s : sub) acc*=s;
    return acc;
}

template<typename T>
std::vector<T> concat(const std::vector<T> & v1, const std::vector<T> & v2){
    auto res = v1;
    res.insert(res.end(), v2.begin(), v2.end());
    return res;
}

template<typename T, auto N, auto M>
auto concat(const std::array<T, N>&a1, const std::array<T, M>&a2){
    std::array<T, N+M> result;
    std::copy (a1.cbegin(), a1.cend(), result.begin());
    std::copy (a2.cbegin(), a2.cend(), result.begin() + N);
    return result;
}

template<typename T, auto N>
auto repeat(const std::array<T, N>&a){
    return concat(a,a);
}

auto group(const std::vector<long> & dim, const std::vector<size_t> &pattern ){
    std::vector<long> res(pattern.size());
    size_t dim_offset = 0;
    for (size_t i = 0; i < pattern.size(); i++){
        long product = 1;
        for(size_t j = 0; j < pattern[i]; j++)
            product *= dim[dim_offset+j];
        res[i] = product;
        dim_offset += pattern[i];
    }
    return res;
}


template<auto N, auto M>
auto group(const std::array<long,N> & dim, const std::array<size_t,M> &pattern ){
    std::array<long,M> res{};
    size_t dim_offset = 0;
    for (size_t i = 0; i < M; i++){
        long product = 1;
        for(size_t j = 0; j < pattern[i]; j++)
            product *= dim[dim_offset+j];
        res[i] = product;
        dim_offset += pattern[i];
    }
    return res;
}


Eigen::Tensor<qm::Scalar,2> contract(const Eigen::Tensor<qm::Scalar,2> & m,
                                     const Eigen::Tensor<qm::Scalar,2> & ud,
                                     const std::array<long,4> & shp_mid4,
                                     const std::array<long,4> & shp_udn4,
                                     const std::array<long,6> & shf6,
                                     const Textra::idxlistpair<1> & idx1,
                                     const Textra::idxlistpair<2> & idx2,
                                     const std::array<long,2> & dim2){
    return m.reshape(shp_mid4)
            .contract(ud.reshape(shp_udn4), idx1)
            .contract(ud.reshape(shp_udn4).conjugate(), idx2)
            .shuffle(shf6).reshape(dim2);
}


Eigen::Tensor<qm::Scalar, 2> qm::Gate::exp_internal(const Eigen::Tensor<Scalar, 2> &op_, Scalar alpha) const {
    /* Note fore flbit:
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
    auto op_map = Textra::TensorMatrixMap(op_);
    if(op_map.isDiagonal() and op_map.imag().isZero() and std::real(alpha) == 0) {
        auto minus_i = std::complex<double>(0, -1);

        auto diag =
            op_map.diagonal()
                .unaryViewExpr([&alpha, &minus_i](const Scalar &h) { return std::exp(minus_i * std::fmod(std::imag(-alpha) * std::real(h), 2.0 * M_PI)); })
                .asDiagonal();
        return Textra::MatrixToTensor(diag);
        //        std::cout << "alpha: " << alpha << std::endl;
        //        std::cout << "op:\n" << op << std::endl;
        //        std::cout << "old:\n" << (alpha * Textra::TensorMatrixMap(op_)).exp() << std::endl;
        //        if(not Textra::TensorMatrixMap(op).isApprox((alpha * Textra::TensorMatrixMap(op_)).exp())){
        //            throw std::runtime_error("Mismatch");
        //        }
    } else {
        return Textra::MatrixToTensor((alpha * Textra::TensorMatrixMap(op_)).exp());
    }
}

qm::Gate::Gate(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> & op_, std::vector<size_t> pos_, std::vector<long> dim_) : pos(std::move(pos_)), dim(std::move(dim_)){
    auto dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies());
    if(dim_prod != op_.rows() or dim_prod != op_.cols()) throw std::logic_error(fmt::format("dim {} not compatible with matrix dimensions {} x {}", dim, op_.rows(), op_.cols()));
    if(pos.size() != dim.size()) throw std::logic_error(fmt::format("pos.size() {} != dim.size() {}", pos, dim));
    op = Textra::MatrixTensorMap(op_);
}


qm::Gate::Gate(const Eigen::Tensor<Scalar, 2> &op_, std::vector<size_t> pos_, std::vector<long>  dim_) : op(op_), pos(std::move(pos_)), dim(std::move(dim_)){
    auto dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies());
    if(dim_prod != op_.dimension(0) or dim_prod != op_.dimension(1)) throw std::logic_error(fmt::format("dim {} not compatible with matrix dimensions {} x {}", dim, op_.dimension(0), op_.dimension(1)));
    if(pos.size() != dim.size()) throw std::logic_error(fmt::format("pos.size() {} != dim.size() {}", pos, dim));
}
qm::Gate::Gate(const Eigen::Tensor<Scalar, 2> &op_, std::vector<size_t> pos_, std::vector<long> dim_, Scalar alpha) : pos(std::move(pos_)),dim(std::move(dim_)){
    auto dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies());
    if(dim_prod != op_.dimension(0) or dim_prod != op_.dimension(1)) throw std::logic_error(fmt::format("dim {} not compatible with matrix dimensions {} x {}", dim, op_.dimension(0), op_.dimension(1)));
    if(pos.size() != dim.size()) throw std::logic_error(fmt::format("pos.size() {} != dim.size() {}", pos, dim));
    op = exp_internal(op_, alpha);
}

void                          qm::Gate::exp_inplace(Scalar alpha) { op = exp_internal(op, alpha); }
qm::Gate                      qm::Gate::exp(Scalar alpha) const { return Gate(op, pos,dim, alpha); }
bool                          qm::Gate::isUnitary(double prec) const { return Textra::TensorMatrixMap(op).isUnitary(prec); }
Eigen::Tensor<qm::Scalar, 2> &qm::Gate::adjoint() const {
    if(adj) return adj.value();
    adj = op.conjugate().shuffle(Textra::array2{1, 0});
    return adj.value();
}

template<auto rank>
std::array<long,rank> qm::Gate::shape() const{
    if constexpr (rank == 1){
        return {op.dimension(0)*op.dimension(1)};
    }
    if constexpr (rank == 2){
        return {op.dimension(0), op.dimension(1)};
    }else{
        if(rank == 2*dim.size()){
            std::array<long,rank> dims{};
            for (size_t i = 0; i < dims.size(); i++) dims[i] = dim[i % dim.size()];
            return dims;
        }else{
            throw std::range_error(fmt::format("Can't compute shape of rank {} for gate with pos {} and dim {}", rank, pos, dim));
        }
    }
}
template std::array<long,1>  qm::Gate::shape() const;
template std::array<long,2>  qm::Gate::shape() const;
template std::array<long,4>  qm::Gate::shape() const;
template std::array<long,6>  qm::Gate::shape() const;
template std::array<long,8>  qm::Gate::shape() const;
template std::array<long,10>  qm::Gate::shape() const;
template std::array<long,12>  qm::Gate::shape() const;
template std::array<long,14>  qm::Gate::shape() const;
template std::array<long,16>  qm::Gate::shape() const;


std::vector<size_t> qm::Gate::idx() const{
    std::vector<size_t> idx(pos.size()*2);
    std::iota(idx.begin(), idx.end(),0);
    return idx;
}

std::vector<size_t> qm::Gate::idx_up() const{
    std::vector<size_t> idx(pos.size());
    std::iota(idx.begin(), idx.end(),0);
    return idx;
}

std::vector<size_t> qm::Gate::idx_dn() const{
    std::vector<size_t> idx(pos.size());
    std::iota(idx.begin(), idx.end(), pos.size());
    return idx;
}

std::vector<size_t> qm::Gate::idx(const std::vector<size_t> &pos_) const{
    std::vector<size_t> idx;
    for(const auto & p_: pos_)
        for(const auto &[i,p]: iter::enumerate(pos) )
            if(p == p_) {
                idx.emplace_back(i);
                idx.emplace_back(i+pos.size());
            }
    return idx;
}


std::vector<size_t> qm::Gate::idx_up(const std::vector<size_t> &pos_) const{
    std::vector<size_t> idx;
    for(const auto & p_: pos_)
        for(const auto &[i,p]: iter::enumerate(pos) )
            if(p == p_) {
                idx.emplace_back(i);
            }
    return idx;
}


std::vector<size_t> qm::Gate::idx_dn(const std::vector<size_t> &pos_) const{
    std::vector<size_t> idx;
    for(const auto & p_: pos_)
        for(const auto &[i,p]: iter::enumerate(pos) )
            if(p == p_) {
                idx.emplace_back(i+pos.size());
            }
    return idx;
}


qm::Gate qm::Gate::insert(const Gate & other) const        { return qm::insert(*this, other);}
qm::Gate qm::Gate::connect_above(const Gate & other) const { return qm::connect(other, *this);}
qm::Gate qm::Gate::connect_under(const Gate & other) const { return qm::connect(*this, other);}

template<auto N>
qm::Gate qm::Gate::trace(const Eigen::array<Eigen::IndexPair<Eigen::Index>, N> & idxpair) const{
    return qm::trace(*this, idxpair);
}
template qm::Gate qm::Gate::trace(const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idxpairs) const;
template qm::Gate qm::Gate::trace(const std::array<Eigen::IndexPair<Eigen::Index>, 2> & idxpairs) const;

qm::Gate qm::Gate::trace_idx(const std::vector<long> & idx_) const{return qm::trace_idx(*this, idx_);}
qm::Gate qm::Gate::trace_pos(const std::vector<size_t> & pos_) const{return qm::trace_pos(*this, pos_);}
qm::Gate qm::Gate::trace_pos(size_t pos_) const{return qm::trace_pos(*this, pos_);}
qm::Scalar qm::Gate::trace() const {return qm::trace(*this);}


qm::Gate qm::insert(const qm::Gate & middle_gate , const qm::Gate & updown_gate){
    std::vector<size_t> pos_isect; // locations that intersect: both on middle and updown gates
    std::vector<size_t> pos_nsect; // locations that do not intersect: not in both on middle and updown gates
    std::set_intersection(middle_gate.pos.begin(),middle_gate.pos.end(),
                          updown_gate.pos.begin(),updown_gate.pos.end(),
                          back_inserter(pos_isect));
    std::set_symmetric_difference(middle_gate.pos.begin(),middle_gate.pos.end(),
                                  updown_gate.pos.begin(),updown_gate.pos.end(),
                                  back_inserter(pos_nsect));
    if(pos_isect.empty()) return middle_gate;
    tools::log->info("Inserting gate pos {} between gates pos {}", middle_gate.pos, updown_gate.pos);
    auto shp_udn4 = updown_gate.shape<4>();
    auto shp_udn2 = updown_gate.shape<2>();
    if(pos_isect.size() == 1 and pos_nsect.size() == 1 and middle_gate.pos.size() == 1 and updown_gate.pos.size() == 2){
        // One common location, one two uncommon. Then this connects a 1-site gate with 2-site gates up and down
        auto shp_mid2 = middle_gate.shape<2>();
        // Decide if this is connects on the left or right leg
        if(middle_gate.pos.front() == pos_isect.front()){
            /*  Right insert
             *
             *      |    |
             *     [  up  ]
             *      |    |                0    1          0
             *      |    |                |    |          |
             *      | [ mid ]      =     [ gate ]  =   [ gate ]
             *      |    |                |    |          |
             *      |    |                2   3           1
             *     [  dn  ]
             *      |    |
             */
            Eigen::Tensor<Scalar,2> op =
                updown_gate.op.reshape(shp_udn4)
                    .contract(middle_gate.op, Textra::idx({3},{0}))
                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({2,3},{0,1}))
                    .reshape(shp_udn2);
            auto pos2 = std::vector<size_t> {updown_gate.pos[0],updown_gate.pos[1]};
            auto dim2 = std::vector<long> {updown_gate.dim[0],updown_gate.dim[1]};
            return qm::Gate{op,pos2,dim2};

        }else{
            /*  Left insert
             *
             *      |    |
             *     [  up  ]
             *      |    |                0    1          0
             *      |    |                |    |          |
             *   [ mid ] |         =     [ gate ]  =   [ gate ]
             *      |    |                |    |          |
             *      |    |                2   3           1
             *     [  dn  ]
             *      |    |
             */
            Eigen::Tensor<Scalar,2> op =
                updown_gate.op.reshape(shp_udn4)
                    .contract(middle_gate.op, Textra::idx({3},{0}))
                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({2,3},{0,1}))
                    .reshape(shp_udn2);
            auto pos2 = std::vector<size_t> {updown_gate.pos[0],updown_gate.pos[1]};
            auto dim2 = std::vector<long> {updown_gate.dim[0],updown_gate.dim[1]};
            return qm::Gate{op,pos2,dim2};
        }

    }
    if(pos_isect.size() == 1 and pos_nsect.size() >= 2 and updown_gate.pos.size() == 2){
        // One common location, so it must be at the edge.
        std::array<long,4> shp_mid4{};
        std::array<long,2> dim2{};
        std::array<Eigen::Index,6> shf6{};
        Textra::idxlistpair<1> idx1;
        Textra::idxlistpair<2> idx2;
        std::vector<size_t> pos;
        std::vector<long> dim;
        size_t merged = pos_nsect.size()-pos_isect.size();
        // Decide if this is connects on the right or left leg
        if(middle_gate.pos.front() == pos_isect.front()) {
            /*  Right insert (Free legs in mid are merged)
             *
             *      |    |     |
             *     [  up  ]    |
             *      |    |     |           0   1   2              0
             *      |    |     |           |   |   |              |
             *      |   [  mid  ]   =    [   gate    ]  =   [   gate    ]
             *      |    |     |           |   |   |              |
             *      |    |     |           3   4   5              1
             *     [  dn  ]    |
             *      |    |     |
             */
            idx1  = Textra::idx({0},{3});
            idx2  = Textra::idx({5,1},{0,1});
            shf6  =  std::array<Eigen::Index,6>{2,3,0,4,5,1};
            pos   = concat(updown_gate.pos, subset(middle_gate.pos, 1, merged));
            dim   = concat(updown_gate.dim, subset(middle_gate.dim, 1, merged));
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
            dim2  = repeat(std::array<long,1>{shp_mid4[1] * shp_udn2[0]});
        }else{
            /*  Left insert (Free legs in mid are merged)
             *
             *     |     |    |
             *     |    [  up  ]
             *     |     |    |           0   1   2              0
             *     |     |    |           |   |   |              |
             *    [  mid  ]   |    =    [   gate    ]  =   [   gate    ]
             *     |     |    |           |   |   |              |
             *     |     |    |           3   4   5              1
             *     |    [  dn  ]
             *     |     |    |
             */
            idx1  = Textra::idx({1},{2});
            idx2  = Textra::idx({2,5},{0,1});
            shf6  =  std::array<Eigen::Index,6>{0,2,3,1,4,5};
            pos   = concat(subset(middle_gate.pos, 0, merged), updown_gate.pos);
            dim   = concat(subset(middle_gate.dim, 0, merged), updown_gate.dim);
            if(pos_nsect.size() == 1) shp_mid4 = group(middle_gate.shape<2>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 2) shp_mid4 = group(middle_gate.shape<4>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 3) shp_mid4 = group(middle_gate.shape<6>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 4) shp_mid4 = group(middle_gate.shape<8>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 5) shp_mid4 = group(middle_gate.shape<10>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 6) shp_mid4 = group(middle_gate.shape<12>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 7) shp_mid4 = group(middle_gate.shape<14>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 8) shp_mid4 = group(middle_gate.shape<16>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 9) shp_mid4 = group(middle_gate.shape<18>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 10) shp_mid4 = group(middle_gate.shape<20>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 11) shp_mid4 = group(middle_gate.shape<22>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 12) shp_mid4 = group(middle_gate.shape<24>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 13) shp_mid4 = group(middle_gate.shape<26>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 14) shp_mid4 = group(middle_gate.shape<28>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 15) shp_mid4 = group(middle_gate.shape<30>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            if(pos_nsect.size() == 16) shp_mid4 = group(middle_gate.shape<32>(), std::array<size_t, 4>{merged, 1, merged, 1 });
            dim2  = repeat(std::array<long,1>{shp_mid4[0] * shp_udn2[0]});
        }
        auto op = contract(middle_gate.op, updown_gate.op, shp_mid4, shp_udn4, shf6, idx1, idx2, dim2);
        return qm::Gate{op,pos,dim};
    }
    if(pos_isect.size() == 2 and pos_nsect.size() == 3 and middle_gate.pos.size() == 5 and updown_gate.pos.size() == 2){
        // All the updown locations get contracted
        auto shp_mid10 = middle_gate.shape<10>();
        // Decide if this is connects on the left, right or somewhere in the center
        if(middle_gate.pos.front() == pos_isect.front()){
            /*  Left insert
             *
             *           |  |  |  |  |
             *          [ up ] |  |  |                0   1   2   3   4              0
             *           |  |  |  |  |                |   |   |   |   |              |
             *          [     mid     ]       =     [       gate       ]   =   [   gate    ]
             *           |  |  |  |  |               |   |   |    |   |              |
             *          [ dn ] |  |  |               5   6   7    8   9              1
             *           |  |  |  |  |
             *
             *     Merge the non-contracted mid legs to simplify the insert contraction
             */
            auto shp_mid4 = group(shp_mid10, std::array<size_t,4>{2,3,2,3});
            auto dim2 = middle_gate.shape<2>();
            Eigen::Tensor<Scalar,2> op =
                updown_gate.op
                    .contract(middle_gate.op.reshape(shp_mid4), Textra::idx({1},{0}))
                    .contract(updown_gate.op.conjugate()      , Textra::idx({2},{0}))
                    .shuffle(Textra::array4{0,1,3,2})
                    .reshape(dim2);
            return qm::Gate{op,middle_gate.pos,middle_gate.dim};
        }else if(middle_gate.pos.back() == pos_isect.back()){
            /*  Right insert
             *
             *           |  |  |  |  |
             *           |  |  | [ up ]               0   1   2   3   4              0
             *           |  |  |  |  |                |   |   |   |   |              |
             *          [     mid     ]       =     [       gate       ]   =   [   gate    ]
             *           |  |  |  |  |               |   |   |    |   |              |
             *           |  |  | [ dn ]              5   6   7    8   9              1
             *           |  |  |  |  |
             *
             *     Merge the non-contracted mid legs to simplify the insert contraction
             */
            auto shp_mid4 = group(shp_mid10, std::array<size_t,4>{3,2,3,2});
            auto dim2 = middle_gate.shape<2>();
            Eigen::Tensor<Scalar,2> op =
                updown_gate.op
                    .contract(middle_gate.op.reshape(shp_mid4), Textra::idx({1},{1}))
                    .contract(updown_gate.op.conjugate(), Textra::idx({3},{0}))
                    .shuffle(Textra::array4{1,0,2,3})
                    .reshape(dim2);
            return qm::Gate{op,middle_gate.pos,middle_gate.dim};
        }
        else if(middle_gate.pos[1] == pos_isect.front()){
            /*  Center insert v1
             *
             *           |  |  |  |  |
             *           | [ up ] |  |                0   1   2   3   4              0
             *           |  |  |  |  |                |   |   |   |   |              |
             *          [     mid     ]       =     [       gate       ]   =   [   gate    ]
             *           |  |  |  |  |               |   |   |    |   |              |
             *           | [ dn ] |  |               5   6   7    8   9              1
             *           |  |  |  |  |
             *
             *
             */
            auto shp_mid6 = group(shp_mid10, std::array<size_t,6>{1,2,2,1,2,2});
            auto dim2 = middle_gate.shape<2>();
            Eigen::Tensor<Scalar,2> op =
                updown_gate.op.reshape(shp_udn2)
                    .contract(middle_gate.op.reshape(shp_mid6), Textra::idx({1},{1}))
                    .contract(updown_gate.op.conjugate(), Textra::idx({4},{0}))
                    .shuffle(Textra::array6{1,0,2,3,5,4})
                    .reshape(dim2);
            return qm::Gate{op,middle_gate.pos,middle_gate.dim};
        }
        else if(middle_gate.pos[2] == pos_isect.front()){
            /*  Center insert v2
             *
             *           |  |  |  |  |
             *           |  | [ up ] |                0   1   2   3   4              0
             *           |  |  |  |  |                |   |   |   |   |              |
             *          [     mid     ]       =     [       gate       ]   =   [   gate    ]
             *           |  |  |  |  |               |   |   |    |   |              |
             *           |  | [ dn ] |               5   6   7    8   9              1
             *           |  |  |  |  |
             *
             *
             */
            auto shp_mid6 = group(shp_mid10, std::array<size_t,6>{2,2,1,2,2,1});
            auto dim2 = middle_gate.shape<2>();
            Eigen::Tensor<Scalar,2> op =
                updown_gate.op.reshape(shp_udn2)
                    .contract(middle_gate.op.reshape(shp_mid6), Textra::idx({1},{1}))
                    .contract(updown_gate.op.conjugate(), Textra::idx({4},{0}))
                    .shuffle(Textra::array6{1,0,2,3,5,4})
                    .reshape(dim2);
            return qm::Gate{op,middle_gate.pos,middle_gate.dim};
        }
        else throw std::runtime_error(fmt::format("Unhandled case: pos_isect {} | pos_nsect {}", pos_isect, pos_nsect));
    }
    if(pos_isect.size() == 2 and pos_nsect.size() == 5 and middle_gate.pos.size() == 7 and updown_gate.pos.size() == 2){
        // All the updown locations get contracted
        auto shp_mid14 = middle_gate.shape<14>();
        size_t usize = updown_gate.pos.size();
        size_t msize = middle_gate.pos.size();
        std::array<long,4> shp_mid4{};
        std::array<long,6> shp_mid6{};
        std::array<long,2> dim2{};
        std::array<Eigen::Index,6> shf6{};
        std::array<Eigen::Index,4> shf4{};
        Textra::idxlistpair<1> idx_up, idx_dn;
        Eigen::Tensor<Scalar,2> op;
        auto offset = static_cast<size_t>(std::distance(middle_gate.pos.begin(), find(middle_gate.pos.begin(), middle_gate.pos.end(), pos_isect.front())));

        // Decide if this is connects on the left, right or somewhere in the center
        if(offset == 0 or offset == 5){
            if(offset == 0){
                /*  Insert at offset 0
                 *
                 *           |  |  |  |  |  |  |
                 *          [ up ] |  |  |  |  |
                 *           |  |  |  |  |  |  |
                 *          [        mid        ]
                 *           |  |  |  |  |  |  |
                 *          [ dn ] |  |  |  |  |
                 *           |  |  |  |  |  |  |
                 *
                 */
                shp_mid4 = group(shp_mid14, std::array<size_t,4>{2,5,2,5});
                dim2 = middle_gate.shape<2>();
                idx_up = Textra::idx({1},{0});
                idx_dn = Textra::idx({2},{0});
                shf4 = Textra::array4{0,1,3,2};
            }else if(offset == 5){
                /*  Insert at offset 5
                 *
                 *           |  |  |  |  |  |  |
                 *           |  |  |  |  | [ up ]
                 *           |  |  |  |  |  |  |
                 *          [        mid        ]
                 *           |  |  |  |  |  |  |
                 *           |  |  |  |  | [ dn ]
                 *           |  |  |  |  |  |  |
                 *
                 */
                shp_mid4 = group(shp_mid14, std::array<size_t,4>{5,2,5,2});
                dim2 = middle_gate.shape<2>();
                idx_up = Textra::idx({1},{1});
                idx_dn = Textra::idx({3},{0});
                shf4 = Textra::array4{1,0,2,3};
            }
            op = updown_gate.op
                .contract(middle_gate.op.reshape(shp_mid4), idx_up)
                .contract(updown_gate.op.conjugate(), idx_dn)
                .shuffle(shf4)
                .reshape(dim2);
        }else{
            /*  Insert at offset 1 to 4
             *
             *           |  |  |  |  |  |  |
             *           | [ up ] |  |  |  |
             *           |  |  |  |  |  |  |
             *          [        mid        ]
             *           |  |  |  |  |  |  |
             *           | [ dn ] |  |  |  |
             *           |  |  |  |  |  |  |
             *
             */
            shp_mid6 = group(shp_mid14, repeat(std::array<size_t,3>{offset,usize,msize-usize-offset}));
            dim2 = middle_gate.shape<2>();
            idx_up = Textra::idx({1},{1});
            idx_dn = Textra::idx({4},{0});
            shf6 = Textra::array6{1,0,2,3,5,4};
            op = updown_gate.op
                .contract(middle_gate.op.reshape(shp_mid6), idx_up)
                .contract(updown_gate.op.conjugate(), idx_dn)
                .shuffle(shf6)
                .reshape(dim2);
        }

        return qm::Gate{op,middle_gate.pos,middle_gate.dim};

    }
    if(pos_isect.size() == 2 and pos_nsect.size() == 7 and middle_gate.pos.size() == 9 and updown_gate.pos.size() == 2){
        // All the updown locations get contracted
        auto shp_mid18 = middle_gate.shape<18>();
        size_t usize = updown_gate.pos.size();
        size_t msize = middle_gate.pos.size();
        std::array<long,4> shp_mid4{};
        std::array<long,6> shp_mid6{};
        std::array<long,2> dim2{};
        std::array<Eigen::Index,6> shf6{};
        std::array<Eigen::Index,4> shf4{};
        Textra::idxlistpair<1> idx_up, idx_dn;
        Eigen::Tensor<Scalar,2> op;
        auto offset = static_cast<size_t>(std::distance(middle_gate.pos.begin(), find(middle_gate.pos.begin(), middle_gate.pos.end(), pos_isect.front())));

        // Decide if this is connects on the left, right or somewhere in the center
        if(offset == 0 or offset == 5){
            if(offset == 0){
                /*  Insert at offset 0
                 *
                 *           |  |  |  |  |  |  |
                 *          [ up ] |  |  |  |  |
                 *           |  |  |  |  |  |  |
                 *          [        mid        ]
                 *           |  |  |  |  |  |  |
                 *          [ dn ] |  |  |  |  |
                 *           |  |  |  |  |  |  |
                 *
                 */
                shp_mid4 = group(shp_mid18, std::array<size_t,4>{2,7,2,7});
                dim2 = middle_gate.shape<2>();
                idx_up = Textra::idx({1},{0});
                idx_dn = Textra::idx({2},{0});
                shf4 = Textra::array4{0,1,3,2};
            }else if(offset == 5){
                /*  Insert at offset 5
                 *
                 *           |  |  |  |  |  |  |
                 *           |  |  |  |  | [ up ]
                 *           |  |  |  |  |  |  |
                 *          [        mid        ]
                 *           |  |  |  |  |  |  |
                 *           |  |  |  |  | [ dn ]
                 *           |  |  |  |  |  |  |
                 *
                 */
                shp_mid4 = group(shp_mid18, std::array<size_t,4>{7,2,7,2});
                dim2 = middle_gate.shape<2>();
                idx_up = Textra::idx({1},{1});
                idx_dn = Textra::idx({3},{0});
                shf4 = Textra::array4{1,0,2,3};
            }
            op = updown_gate.op
                .contract(middle_gate.op.reshape(shp_mid4), idx_up)
                .contract(updown_gate.op.conjugate(), idx_dn)
                .shuffle(shf4)
                .reshape(dim2);
        }else{
            /*  Insert at offset 1 to 7
             *
             *           |  |  |  |  |  |  |
             *           | [ up ] |  |  |  |
             *           |  |  |  |  |  |  |
             *          [        mid        ]
             *           |  |  |  |  |  |  |
             *           | [ dn ] |  |  |  |
             *           |  |  |  |  |  |  |
             *
             */
            shp_mid6 = group(shp_mid18, repeat(std::array<size_t,3>{offset,usize,msize-usize-offset}));
            dim2 = middle_gate.shape<2>();
            idx_up = Textra::idx({1},{1});
            idx_dn = Textra::idx({4},{0});
            shf6 = Textra::array6{1,0,2,3,5,4};
            op = updown_gate.op
                .contract(middle_gate.op.reshape(shp_mid6), idx_up)
                .contract(updown_gate.op.conjugate(), idx_dn)
                .shuffle(shf6)
                .reshape(dim2);
        }

        return qm::Gate{op,middle_gate.pos,middle_gate.dim};

    }
    throw std::runtime_error(fmt::format("Case not implemented: pos_isect {} | pos_nsect {}", pos_isect, pos_nsect));
}

qm::Gate qm::connect(const qm::Gate & dn_gate , const qm::Gate & up_gate){
    std::vector<size_t> pos_isect; // locations that intersect: both on middle and updown gates
    std::vector<size_t> pos_nsect; // locations that do not intersect: not in both on middle and updown gates
    std::set_intersection(dn_gate.pos.begin(),dn_gate.pos.end(),
                          up_gate.pos.begin(),up_gate.pos.end(),
                          back_inserter(pos_isect));
    std::set_symmetric_difference(dn_gate.pos.begin(),dn_gate.pos.end(),
                                  up_gate.pos.begin(),up_gate.pos.end(),
                                  back_inserter(pos_nsect));
    if(pos_isect.empty()) return dn_gate;
    if(pos_isect.size() == 1 and pos_nsect.size() == 2){
        // One common location, and two uncommon. Then this connects two 2-site gates
        auto dim_dn4 = dn_gate.shape<4>();
        auto dim_dn2 = dn_gate.shape<2>();
        auto dim_up4 = up_gate.shape<4>();
        auto dim_up2 = up_gate.shape<2>();
        // Decide if this is connects on the right or right leg
        bool right = dn_gate.pos.front() == pos_isect.front();
        if(right){
            /*  Right connection
             *
             *      |    |     |           0   1   2              0
             *     [  up  ]    |           |   |   |              |
             *      |    |     |   =     [   gate    ]  =   [   gate    ]
             *      |   [  dn  ]           |   |   |              |
             *      |    |     |           3   4   5              1
             *
             */

            auto dim2 = std::array<long,2> {dim_up2[0] * dim_dn4[1], dim_up2[1] * dim_dn4[2]};
            Eigen::Tensor<Scalar,2> op =
                dn_gate.op.reshape(dim_dn4)
                    .contract(up_gate.op.reshape(dim_up4), Textra::idx({0},{3}))
                    .shuffle(Textra::array6{3,4,0,5,1,2})
                    .reshape(dim2);

            auto pos3 = std::vector<size_t> {up_gate.pos[0],up_gate.pos[1], dn_gate.pos[1]};
            auto dim3 = std::vector<long> {up_gate.dim[0],up_gate.dim[1], dn_gate.dim[1]};
            return qm::Gate{op,pos3,dim3};
        }else{
            /*  Left connection
             *
             *     |     |    |             0   1   2              0
             *     |    [  up  ]            |   |   |              |
             *     |     |    |      =    [   gate    ]  =   [   gate    ]
             *    [  dn  ]    |             |   |   |              |
             *     |     |    |             3   4   5              1
             */
            auto dim2 = std::array<long,2> {dim_dn4[0]* dim_up2[0], dim_dn4[2] * dim_up2[1] };
            Eigen::Tensor<Scalar,2> op =
                dn_gate.op.reshape(dim_dn4)
                    .contract(up_gate.op.reshape(dim_up4), Textra::idx({1},{2}))
                    .shuffle(Textra::array6{0,3,4,1,2,5})
                    .reshape(dim2);

            auto pos3 = std::vector<size_t> {dn_gate.pos[0],up_gate.pos[0], up_gate.pos[1]};
            auto dim3 = std::vector<long> {dn_gate.dim[0],up_gate.dim[0], up_gate.dim[1]};
            return qm::Gate{op,pos3,dim3};
        }
    }
    throw std::runtime_error(fmt::format("Case not implemented: pos_isect {} | pos_nsect {}", pos_isect, pos_nsect));
}


qm::Gate qm::trace_idx(const qm::Gate & gate , const std::vector<long> & idx){
    if(idx.size() == 2) return qm::trace(gate, Textra::idx({idx[0]},{idx[1]}));
    if(idx.size() == 4) return qm::trace(gate, Textra::idx({idx[0],idx[1]},{idx[2],idx[3]}));
    throw std::runtime_error(fmt::format("Tracing {} indices is not implemented", idx.size()));
}

qm::Gate qm::trace_pos(const qm::Gate & gate , const std::vector<size_t> & pos){
    if(pos.size() <= 2){
        auto idx = gate.idx(pos);
        std::vector<long> idx_long(idx.begin(),idx.end());
        return qm::trace_idx(gate, idx_long);
    }
    qm::Gate tmp = gate;
    for(auto &p : pos){
        tmp = tmp.trace_pos(std::vector<size_t>{p});
    }
    return tmp;
}

qm::Gate qm::trace_pos(const qm::Gate & gate ,size_t pos){
    return qm::trace_pos(gate, gate.idx(std::vector<size_t>{pos}));
}
qm::Scalar qm::trace(const qm::Gate & gate)  {
    qm::Gate t = qm::trace_pos(gate, gate.pos);
    return t.op(0);
}

template<auto N>
qm::Gate qm::trace(const qm::Gate & gate ,  const std::array<Eigen::IndexPair<Eigen::Index>, N>  & idxpairs){
    // Compute the remaining indices
    auto idx = gate.idx();
    for(auto & pair : idxpairs){
        idx.erase(std::remove(idx.begin(), idx.end(), pair.first), idx.end());
        idx.erase(std::remove(idx.begin(), idx.end(), pair.second), idx.end());
    }

    // Extract the remaining positions and dimensions
    std::vector<size_t> pos;
    std::vector<long> dim;
    for(const auto & i : idx) {
        pos.emplace_back(gate.pos[i]);
        dim.emplace_back(gate.dim[i]);
    }

    // Assuming the gate is symmetric, we can compute the rank2 dimensions for the storage of the gate op
    long dim_prod = std::accumulate(std::begin(dim), std::end(dim), 1, std::multiplies()); // Product of all dimensions of the remaining top legs of the gate
    std::array<long,2> dim2{dim_prod, dim_prod};

    using T2 = Eigen::Tensor<Scalar,2>;
    using T4 = Eigen::Tensor<Scalar,4>;
    using T6 = Eigen::Tensor<Scalar,6>;
    using T8 = Eigen::Tensor<Scalar,8>;
    using T10 = Eigen::Tensor<Scalar,10>;
    using T12 = Eigen::Tensor<Scalar,12>;
    using T14 = Eigen::Tensor<Scalar,14>;
    using T16 = Eigen::Tensor<Scalar,16>;
    using T18 = Eigen::Tensor<Scalar,18>;
    using T20 = Eigen::Tensor<Scalar,20>;
    using T22 = Eigen::Tensor<Scalar,22>;
    using T24 = Eigen::Tensor<Scalar,24>;
    using T26 = Eigen::Tensor<Scalar,26>;
    Eigen::Tensor<Scalar,2> op_traced;
    // Trace
    if constexpr (N == 1)
        if(gate.dim.size() == 1){
            auto shape = gate.shape<2>();
            op_traced = Textra::trace(static_cast<T2>(gate.op.reshape(shape)), idxpairs).reshape(dim2);
            return Gate{op_traced,pos,dim};
        }
    if(gate.dim.size() == 2) op_traced = Textra::trace(static_cast<T4>(gate.op.reshape(gate.shape<4>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 3) op_traced = Textra::trace(static_cast<T6>(gate.op.reshape(gate.shape<6>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 4) op_traced = Textra::trace(static_cast<T8>(gate.op.reshape(gate.shape<8>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 5) op_traced = Textra::trace(static_cast<T10>(gate.op.reshape(gate.shape<10>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 6) op_traced = Textra::trace(static_cast<T12>(gate.op.reshape(gate.shape<12>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 7) op_traced = Textra::trace(static_cast<T14>(gate.op.reshape(gate.shape<14>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 8) op_traced = Textra::trace(static_cast<T16>(gate.op.reshape(gate.shape<16>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 9) op_traced = Textra::trace(static_cast<T18>(gate.op.reshape(gate.shape<18>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 10) op_traced = Textra::trace(static_cast<T20>(gate.op.reshape(gate.shape<20>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 11) op_traced = Textra::trace(static_cast<T22>(gate.op.reshape(gate.shape<22>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 12) op_traced = Textra::trace(static_cast<T24>(gate.op.reshape(gate.shape<24>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 13) op_traced = Textra::trace(static_cast<T24>(gate.op.reshape(gate.shape<26>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 14) op_traced = Textra::trace(static_cast<T24>(gate.op.reshape(gate.shape<28>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 15) op_traced = Textra::trace(static_cast<T24>(gate.op.reshape(gate.shape<30>())), idxpairs).reshape(dim2);
    if(gate.dim.size() == 16) op_traced = Textra::trace(static_cast<T24>(gate.op.reshape(gate.shape<32>())), idxpairs).reshape(dim2);
    return Gate{op_traced,pos,dim};
    throw std::runtime_error(fmt::format("Trace not implemented: N == {} | dim.size() == {}", N, gate.dim.size()));
}


template qm::Gate qm::trace(const qm::Gate & gate , const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idxpairs);
template qm::Gate qm::trace(const qm::Gate & gate , const std::array<Eigen::IndexPair<Eigen::Index>, 2> & idxpairs);

//
//using T2 = Eigen::Tensor<Scalar,2>;
//using T3 = Eigen::Tensor<Scalar,3>;
//using T4 = Eigen::Tensor<Scalar,4>;
//using T5 = Eigen::Tensor<Scalar,5>;
//using T6 = Eigen::Tensor<Scalar,6>;
//using T7 = Eigen::Tensor<Scalar,7>;
//using T8 = Eigen::Tensor<Scalar,8>;
//using T9 = Eigen::Tensor<Scalar,9>;
//using T10 = Eigen::Tensor<Scalar,10>;
//using T11 = Eigen::Tensor<Scalar,11>;
//using T12 = Eigen::Tensor<Scalar,12>;
//Eigen::Tensor<Scalar,2> op_traced;
//// Trace
//if constexpr (N == 1)
//if(gate.dim.size() == 1) op_traced = Textra::trace(static_cast<T2>(gate.op.reshape(gate.shape<2>())), idxpairs).reshape(dim2);
//
//if(gate.dim.size() == 2) op_traced = Textra::trace(static_cast<T2>(gate.op.reshape(gate.shape<2>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 3) op_traced = Textra::trace(static_cast<T3>(gate.op.reshape(gate.shape<3>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 4) op_traced = Textra::trace(static_cast<T4>(gate.op.reshape(gate.shape<4>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 5) op_traced = Textra::trace(static_cast<T5>(gate.op.reshape(gate.shape<5>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 6) op_traced = Textra::trace(static_cast<T6>(gate.op.reshape(gate.shape<6>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 7) op_traced = Textra::trace(static_cast<T7>(gate.op.reshape(gate.shape<7>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 8) op_traced = Textra::trace(static_cast<T8>(gate.op.reshape(gate.shape<8>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 9) op_traced = Textra::trace(static_cast<T9>(gate.op.reshape(gate.shape<9>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 10) op_traced = Textra::trace(static_cast<T10>(gate.op.reshape(gate.shape<10>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 11) op_traced = Textra::trace(static_cast<T11>(gate.op.reshape(gate.shape<11>())), idxpairs).reshape(dim2);
//if(gate.dim.size() == 12) op_traced = Textra::trace(static_cast<T12>(gate.op.reshape(gate.shape<12>())), idxpairs).reshape(dim2);



//    if(pos_isect.size() == 1 and pos_nsect.size() == 2 and middle_gate.pos.size() == 2 and updown_gate.pos.size() == 2){
//        // One common location, and two uncommon. Then this connects two 2-site gates
//        auto shp_mid4 = middle_gate.shape<4>();
//        // Decide if this is connects on the right or right leg
//        if(middle_gate.pos.front() == pos_isect.front()){
//            /*  Right insert
//             *
//             *      |    |     |
//             *     [  up  ]    |
//             *      |    |     |                0   1   2              0
//             *      |    |     |                |   |   |              |
//             *      |   [  mid  ]       =     [   gate    ]  =   [   gate    ]
//             *      |    |     |                |   |   |              |
//             *      |    |     |                3   4   5              1
//             *     [  dn  ]    |
//             *      |    |     |
//             */
//            auto dim2 = std::array<long,2> {shp_udn2[0] * shp_mid4[1], shp_udn2[1] * shp_mid4[3]};
//            Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({0},{3}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({5,1},{0,1}))
//                    .shuffle(Textra::array6{2,3,0,4,5,1})
//                    .reshape(dim2);
//
//            auto pos3 = std::vector<size_t> {updown_gate.pos[0],updown_gate.pos[1], middle_gate.pos[1]};
//            auto dim3 = std::vector<long> {updown_gate.dim[0],updown_gate.dim[1], middle_gate.dim[1]};
//            return qm::Gate{op,pos3,dim3};
//        }else{
//            /*  Left insert
//             *
//             *     |     |    |
//             *     |    [  up  ]
//             *     |     |    |           0   1   2              0
//             *     |     |    |           |   |   |              |
//             *    [  mid  ]   |    =    [   gate    ]  =   [   gate    ]
//             *     |     |    |           |   |   |              |
//             *     |     |    |           3   4   5              1
//             *     |    [  dn  ]
//             *     |     |    |
//             */
//            auto dim2 = std::array<long,2> {shp_mid4[0]* shp_udn2[0], shp_mid4[2] * shp_udn2[1] };
//            Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({1},{2}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({2,5},{0,1}))
//                    .shuffle(Textra::array6{0,2,3,1,4,5})
//                    .reshape(dim2);
//
//            auto pos3 = std::vector<size_t> {middle_gate.pos[0],updown_gate.pos[0], updown_gate.pos[1]};
//            auto dim3 = std::vector<long> {middle_gate.dim[0],updown_gate.dim[0], updown_gate.dim[1]};
//            return qm::Gate{op,pos3,dim3};
//        }
//    }
//    if(pos_isect.size() == 1 and pos_nsect.size() == 3 and middle_gate.pos.size() == 3 and updown_gate.pos.size() == 2){
//        // One common location, and three uncommon. Then this connects a 3-site gate between 2-site gates
//        auto shp_mid6 = middle_gate.shape<6>();
//        // Decide if this is connects on the right or right leg
//        if(middle_gate.pos.front() == pos_isect.front()){
//            /*  Right insert
//             *
//             *      |    |    |    |
//             *     [  up  ]   |    |
//             *      |    |    |    |                0   1   2   3              0
//             *      |    |    |    |                |   |   |   |              |
//             *      |   [    mid    ]       =     [    gate      ]  =   [   gate    ]
//             *      |    |    |    |               |   |   |    |              |
//             *      |    |    |    |               4   5   6    7              1
//             *     [  dn  ]   |    |
//             *      |    |    |    |
//             *
//             *      Merge the free legs of mid to simplify the contraction! (See shp_mid4)
//             *      Doing so allows us to use the same contraction used in the 2-2 legged left insert
//             *
//             */
//            auto shp_mid4 = group(shp_mid6, std::array<size_t,4>{1,2,1,2});
//            auto dim2 = std::array<long,2> {shp_udn2[0] * shp_mid4[1], shp_udn2[1] * shp_mid4[3]};
//
//            Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({0},{3}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({5,1},{0,1}))
//                    .shuffle(Textra::array6{2,3,0,4,5,1})
//                    .reshape(dim2);
//
//            auto pos4 = std::vector<size_t> {updown_gate.pos[0],updown_gate.pos[1], middle_gate.pos[1], middle_gate.pos[2]};
//            auto dim4 = std::vector<long> {updown_gate.dim[0],updown_gate.dim[1], middle_gate.dim[1], middle_gate.dim[2]};
//            return qm::Gate{op,pos4,dim4};
//        }else{
//            /*  Left insert
//             *
//             *     |    |    |    |
//             *     |    |   [  up  ]
//             *     |    |    |    |           0   1   2   3              0
//             *     |    |    |    |           |   |   |   |              |
//             *    [    mid    ]   |    =    [    gate      ]  =   [   gate    ]
//             *     |     |   |    |          |   |   |    |              |
//             *     |     |   |    |          4   5   6    7              1
//             *     |     |  [  dn  ]
//             *     |     |   |    |
//             *
//             *      Merge the free legs of mid to simplify the contraction! (See shp_mid4)
//             *      Doing so allows us to use the same contraction used in the 2-2 legged left insert
//             *
//             */
//            auto shp_mid4 = group(shp_mid6, std::array<size_t,4>{2,1,2,1});
//            auto dim2 = std::array<long,2> {shp_mid4[0] * shp_udn2[0], shp_mid4[2] * shp_udn2[1] };
//             Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({1},{2}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({2,5},{0,1}))
//                    .shuffle(Textra::array6{0,2,3,1,4,5})
//                    .reshape(dim2);
//
//            auto pos4 = std::vector<size_t> {middle_gate.pos[0], middle_gate.pos[1],updown_gate.pos[0], updown_gate.pos[1]};
//            auto dim4 = std::vector<long> {middle_gate.dim[0], middle_gate.dim[1],updown_gate.dim[0], updown_gate.dim[1]};
//            return qm::Gate{op,pos4,dim4};
//        }
//    }
//    if(pos_isect.size() == 1 and pos_nsect.size() == 4 and middle_gate.pos.size() == 4 and updown_gate.pos.size() == 2){
//        // One common location, and 4 uncommon. Then this connects a 4-site gate between 2-site gates
//        auto shp_mid8 = middle_gate.shape<8>();
//        // Decide if this is connects on the right or right leg
//        if(middle_gate.pos.front() == pos_isect.front()){
//            /*  Right insert
//             *
//             *      |    |  |   |  |
//             *     [  up  ] |   |  |
//             *      |    |  |   |  |                0   1   2   3   4              0
//             *      |    |  |   |  |                |   |   |   |   |              |
//             *      |   [    mid    ]       =     [       gate       ]   =   [   gate    ]
//             *      |    |  |   |  |               |   |   |    |   |              |
//             *      |    |  |   |  |               5   6   7    8   9              1
//             *     [  dn  ] |   |  |
//             *      |    |  |   |  |
//             *
//             *      Merge the free legs of mid to simplify the contraction! (See shp_mid4)
//             *      Doing so allows us to use the same contraction used in the 2-2 legged left insert
//             *
//             */
//            auto shp_mid4 = group(shp_mid8, std::array<size_t,4>{1,3,1,3});
//            auto dim2 = std::array<long,2> {shp_udn2[0] * shp_mid4[1], shp_udn2[1] * shp_mid4[3]};
//            Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({0},{3}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({5,1},{0,1}))
//                    .shuffle(Textra::array6{2,3,0,4,5,1})
//                    .reshape(dim2);
//            auto pos5 = concat(updown_gate.pos, subset(middle_gate.pos, 1, 3));
//            auto dim5 = concat(updown_gate.dim, subset(middle_gate.dim, 1, 3));
//            return qm::Gate{op,pos5,dim5};
//        }else{
//            /*  Left insert
//             *
//             *     |   |   |   |    |
//             *     |   |   |  [  up  ]
//             *     |   |   |   |    |           0   1   2   3   4              0
//             *     |   |   |   |    |           |   |   |   |   |              |
//             *    [     mid     ]   |    =    [       gate       ]   =   [   gate    ]
//             *     |   |   |   |    |          |   |   |    |   |              |
//             *     |   |   |   |    |          5   6   7    8   9              1
//             *     |   |   |  [  dn  ]
//             *     |   |   |   |    |
//             *
//             *      Merge the free legs of mid to simplify the contraction! (See shp_mid4)
//             *      Doing so allows us to use the same contraction used in the 2-2 legged left insert
//             *
//             */
//            auto shp_mid4 = group(shp_mid8, std::array<size_t,4>{3,1,3,1});
//            auto dim2 = std::array<long,2> {shp_mid4[0] * shp_udn2[0], shp_mid4[2] * shp_udn2[1] };
//            Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({1},{2}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({2,5},{0,1}))
//                    .shuffle(Textra::array6{0,2,3,1,4,5})
//                    .reshape(dim2);
//            auto pos5 = concat(subset(middle_gate.pos, 0, 3), updown_gate.pos);
//            auto dim5 = concat(subset(middle_gate.dim, 0, 3), updown_gate.dim);
//            return qm::Gate{op,pos5,dim5};
//        }
//    }
//    if(pos_isect.size() == 1 and pos_nsect.size() == 5 and middle_gate.pos.size() == 5 and updown_gate.pos.size() == 2){
//        // One common location, and 4 uncommon. Then this connects a 4-site gate between 2-site gates
//        auto shp_mid10 = middle_gate.shape<10>();
//        // Decide if this is connects on the right or right leg
//        bool right = middle_gate.pos.front() == pos_isect.front();
//        if(right){
//            /*  Right insert
//             *
//             *      |    |  |  |  |  |
//             *     [  up  ] |  |  |  |
//             *      |    |  |  |  |  |                0   1   2   3   4   5              0
//             *      |    |  |  |  |  |                |   |   |   |   |   |              |
//             *      |   [     mid     ]       =     [        gate          ]   =   [   gate    ]
//             *      |    |  |  |  |  |               |   |   |   |    |   |              |
//             *      |    |  |  |  |  |               6   7   8   9   10  11              1
//             *     [  dn  ] |  |  |  |
//             *      |    |  |  |  |  |
//             *
//             *      Merge the free legs of mid to simplify the contraction! (See shp_mid4)
//             *      Doing so allows us to use the same contraction used in the 2-2 legged left insert
//             *
//             */
//            auto shp_mid4 = group(shp_mid10, std::array<size_t,4>{1,4,1,4});
//            auto dim2 = std::array<long,2> {shp_udn2[0] * shp_mid4[1], shp_udn2[1] * shp_mid4[3]};
//            Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({0},{3}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({5,1},{0,1}))
//                    .shuffle(Textra::array6{2,3,0,4,5,1})
//                    .reshape(dim2);
//            auto pos6 = concat(updown_gate.pos, subset(middle_gate.pos, 1, 4));
//            auto dim6 = concat(updown_gate.dim, subset(middle_gate.dim, 1, 4));
//            return qm::Gate{op,pos6,dim6};
//        }else{
//            /*  Left insert
//             *
//             *     |  |  |  |   |    |
//             *     |  |  |  |  [  up  ]
//             *     |  |  |  |   |    |            0   1   2   3   4   5              0
//             *     |  |  |  |   |    |            |   |   |   |   |   |              |
//             *    [      mid      ]  |    =     [        gate          ]   =   [   gate    ]
//             *     |  |  |  |   |    |           |   |   |   |    |   |              |
//             *     |  |  |  |   |    |           6   7   8   9   10  11              1
//             *     |  |  |  |  [  dn  ]
//             *     |  |  |  |   |    |
//             *
//             *      Merge the free legs of mid to simplify the contraction! (See shp_mid4)
//             *      Doing so allows us to use the same contraction used in the 2-2 legged left insert
//             *
//             */
//            auto shp_mid4 = group(shp_mid10, std::array<size_t,4>{4,1,4,1});
//            auto dim2 = std::array<long,2> {shp_mid4[0] * shp_udn2[0], shp_mid4[2] * shp_udn2[1] };
//            Eigen::Tensor<Scalar,2> op =
//                middle_gate.op.reshape(shp_mid4)
//                    .contract(updown_gate.op.reshape(shp_udn4), Textra::idx({1},{2}))
//                    .contract(updown_gate.op.reshape(shp_udn4).conjugate(), Textra::idx({2,5},{0,1}))
//                    .shuffle(Textra::array6{0,2,3,1,4,5})
//                    .reshape(dim2);
//            auto pos6 = concat(subset(middle_gate.pos, 0, 4), updown_gate.pos);
//            auto dim6 = concat(subset(middle_gate.dim, 0, 4), updown_gate.dim);
//            return qm::Gate{op,pos6,dim6};
//        }
//    }
//    if(pos_isect.size() == 1 and pos_nsect.size() == 6 and middle_gate.pos.size() == 6 and updown_gate.pos.size() == 2){
//        // One common location, and 4 uncommon. Then this connects a 4-site gate between 2-site gates
//        auto shp_mid12 = middle_gate.shape<12>();
//        std::array<long,4> shp_mid4{};
//        std::array<long,2> dim2{};
//        std::array<Eigen::Index,6> shf6{};
//        Textra::idxlistpair<1> idx1;
//        Textra::idxlistpair<2> idx2;
//        std::vector<size_t> pos;
//        std::vector<long> dim;
//        // Decide if this is connects on the right or left leg
//        if(middle_gate.pos.front() == pos_isect.front()) {
//            shp_mid4 = group(shp_mid12, std::array<size_t, 4>{1, merged, 1, merged});
//            dim2  = repeat(std::array<long,1>{shp_mid4[1] * shp_udn2[0]});
//            idx1  = Textra::idx({0},{3});
//            idx2  = Textra::idx({5,1},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{2,3,0,4,5,1};
//            pos   = concat(updown_gate.pos, subset(middle_gate.pos, 1, merged));
//            dim   = concat(updown_gate.dim, subset(middle_gate.dim, 1, merged));
//        }else {
//            shp_mid4 = group(shp_mid12, std::array<size_t, 4>{merged, 1, merged, 1});
//            dim2  = repeat(std::array<long,1>{shp_mid4[0] * shp_udn2[0]});
//            idx1  = Textra::idx({1},{2});
//            idx2  = Textra::idx({2,5},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{0,2,3,1,4,5};
//            pos   = concat(subset(middle_gate.pos, 0, merged), updown_gate.pos);
//            dim   = concat(subset(middle_gate.dim, 0, merged), updown_gate.dim);
//        }
//        Eigen::Tensor<Scalar,2> op =
//            middle_gate.op.reshape(shp_mid4)
//                .contract(updown_gate.op.reshape(shp_udn4), idx1)
//                .contract(updown_gate.op.reshape(shp_udn4).conjugate(), idx2)
//                .shuffle(shf6).reshape(dim2);
//        return qm::Gate{op,pos,dim};
//    }
//    if(pos_isect.size() == 1 and pos_nsect.size() == 7 and middle_gate.pos.size() == 7 and updown_gate.pos.size() == 2){
//        // One common location, and 4 uncommon. Then this connects a 4-site gate between 2-site gates
//        auto shp_mid14 = middle_gate.shape<14>();
//        std::array<long,4> shp_mid4{};
//        std::array<long,2> dim2{};
//        std::array<Eigen::Index,6> shf6{};
//        Textra::idxlistpair<1> idx1;
//        Textra::idxlistpair<2> idx2;
//        std::vector<size_t> pos;
//        std::vector<long> dim;
//        // Decide if this is connects on the right or left leg
//        if(middle_gate.pos.front() == pos_isect.front()) {
//            shp_mid4 = group(shp_mid14, std::array<size_t, 4>{1, merged, 1, merged});
//            dim2  = repeat(std::array<long,1>{shp_mid4[1] * shp_udn2[0]});
//            idx1  = Textra::idx({0},{3});
//            idx2  = Textra::idx({5,1},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{2,3,0,4,5,1};
//            pos   = concat(updown_gate.pos, subset(middle_gate.pos, 1, merged));
//            dim   = concat(updown_gate.dim, subset(middle_gate.dim, 1, merged));
//        }else {
//            shp_mid4 = group(shp_mid14, std::array<size_t, 4>{merged, 1, merged, 1});
//            dim2  = repeat(std::array<long,1>{shp_mid4[0] * shp_udn2[0]});
//            idx1  = Textra::idx({1},{2});
//            idx2  = Textra::idx({2,5},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{0,2,3,1,4,5};
//            pos   = concat(subset(middle_gate.pos, 0, merged), updown_gate.pos);
//            dim   = concat(subset(middle_gate.dim, 0, merged), updown_gate.dim);
//        }
//        Eigen::Tensor<Scalar,2> op =
//            middle_gate.op.reshape(shp_mid4)
//                .contract(updown_gate.op.reshape(shp_udn4), idx1)
//                .contract(updown_gate.op.reshape(shp_udn4).conjugate(), idx2)
//                .shuffle(shf6).reshape(dim2);
//        return qm::Gate{op,pos,dim};
//    }
//    if(pos_isect.size() == 1 and pos_nsect.size() == 8 and middle_gate.pos.size() == 8 and updown_gate.pos.size() == 2){
//        // One common location, and 4 uncommon. Then this connects a 4-site gate between 2-site gates
//        auto shp_mid16 = middle_gate.shape<16>();
//        std::array<long,4> shp_mid4{};
//        std::array<long,2> dim2{};
//        std::array<Eigen::Index,6> shf6{};
//        Textra::idxlistpair<1> idx1;
//        Textra::idxlistpair<2> idx2;
//        std::vector<size_t> pos;
//        std::vector<long> dim;
//        // Decide if this is connects on the right or left leg
//        if(middle_gate.pos.front() == pos_isect.front()) {
//            shp_mid4 = group(shp_mid16, std::array<size_t, 4>{1, merged, 1, merged});
//            dim2  = repeat(std::array<long,1>{shp_mid4[1] * shp_udn2[0]});
//            idx1  = Textra::idx({0},{3});
//            idx2  = Textra::idx({5,1},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{2,3,0,4,5,1};
//            pos   = concat(updown_gate.pos, subset(middle_gate.pos, 1, merged));
//            dim   = concat(updown_gate.dim, subset(middle_gate.dim, 1, merged));
//        }else {
//            shp_mid4 = group(shp_mid16, std::array<size_t, 4>{merged, 1, merged, 1});
//            dim2  = repeat(std::array<long,1>{shp_mid4[0] * shp_udn2[0]});
//            idx1  = Textra::idx({1},{2});
//            idx2  = Textra::idx({2,5},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{0,2,3,1,4,5};
//            pos   = concat(subset(middle_gate.pos, 0, merged), updown_gate.pos);
//            dim   = concat(subset(middle_gate.dim, 0, merged), updown_gate.dim);
//        }
//        auto op = contract(middle_gate.op, updown_gate.op, shp_mid4, shp_udn4, shf6, idx1, idx2, dim2);
//        return qm::Gate{op,pos,dim};
//    }
//    if(pos_isect.size() == 1 and pos_nsect.size() == 9 and middle_gate.pos.size() == 9 and updown_gate.pos.size() == 2){
//        // One common location, and 4 uncommon. Then this connects a 4-site gate between 2-site gates
//        auto shp_mid18 = middle_gate.shape<18>();
//        std::array<long,4> shp_mid4{};
//        std::array<long,2> dim2{};
//        std::array<Eigen::Index,6> shf6{};
//        Textra::idxlistpair<1> idx1;
//        Textra::idxlistpair<2> idx2;
//        std::vector<size_t> pos;
//        std::vector<long> dim;
//        // Decide if this is connects on the right or left leg
//        if(middle_gate.pos.front() == pos_isect.front()) {
//            shp_mid4 = group(shp_mid18, std::array<size_t, 4>{1, merged, 1, merged});
//            dim2  = repeat(std::array<long,1>{shp_mid4[1] * shp_udn2[0]});
//            idx1  = Textra::idx({0},{3});
//            idx2  = Textra::idx({5,1},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{2,3,0,4,5,1};
//            pos   = concat(updown_gate.pos, subset(middle_gate.pos, 1, merged));
//            dim   = concat(updown_gate.dim, subset(middle_gate.dim, 1, merged));
//        }else {
//            shp_mid4 = group(shp_mid18, std::array<size_t, 4>{merged, 1, merged, 1});
//            dim2  = repeat(std::array<long,1>{shp_mid4[0] * shp_udn2[0]});
//            idx1  = Textra::idx({1},{2});
//            idx2  = Textra::idx({2,5},{0,1});
//            shf6  =  std::array<Eigen::Index,6>{0,2,3,1,4,5};
//            pos   = concat(subset(middle_gate.pos, 0, merged), updown_gate.pos);
//            dim   = concat(subset(middle_gate.dim, 0, merged), updown_gate.dim);
//        }
//        auto op = contract(middle_gate.op, updown_gate.op, shp_mid4, shp_udn4, shf6, idx1, idx2, dim2);
//        return qm::Gate{op,pos,dim};
//    }
