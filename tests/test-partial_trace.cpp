#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <io/spdlog.h>
#include <math/linalg.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <iostream>

namespace mat{

}
template<typename T>
std::string to_string(const Eigen::EigenBase<T> & m){
//    for(size_t i = 0; i < static_cast<size_t>(m.size()); i++){
//        fmt::print("{}: {}\n",i, m.derived().data()[i]);
//    }
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "  [", "]");
    std::stringstream ss;
    ss << m.derived().format(CleanFmt);;
    return ss.str();
}

template<typename T>
std::string to_string(const Eigen::TensorBase<T,Eigen::ReadOnlyAccessors> & expr){
    typedef Eigen::TensorEvaluator<const Eigen::TensorForcedEvalOp<const T>, Eigen::DefaultDevice> Evaluator;
    typedef typename Evaluator::Dimensions Dimensions;

    // Evaluate the expression if needed
    Eigen::TensorForcedEvalOp<const T> eval = expr.eval();
    Evaluator tensor(eval, Eigen::DefaultDevice());
    tensor.evalSubExprsIfNeeded(NULL);
    Eigen::Index total_size = Eigen::internal::array_prod(tensor.dimensions());

//    for(size_t i = 0; i < static_cast<size_t>(total_size); i++){
//        fmt::print("{}: {}\n",i, tensor.data()[i]);
//    }

    if(total_size > 0){
        Eigen::Index first_dim = tensor.dimensions()[0];
        Eigen::Index other_dim = total_size/first_dim;
        if constexpr(T::NumDimensions == 4) first_dim = tensor.dimensions()[0] * tensor.dimensions()[1];
        if constexpr(T::NumDimensions == 6) first_dim = tensor.dimensions()[0] * tensor.dimensions()[1] *tensor.dimensions()[2];

        other_dim = total_size/first_dim;

        Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "  [", "]");
        std::stringstream ss;
        ss << Eigen::Map<Eigen::Matrix<typename T::Scalar,Eigen::Dynamic,Eigen::Dynamic,Evaluator::Layout>>(tensor.data(),first_dim,other_dim).format(CleanFmt);
        return ss.str();
    }

}


template<typename A, typename B>
bool equal (const A & a, const B & b, double prec = 1e-8){
    if(a.size() != b.size()) return false;
    for(size_t i = 0; i < static_cast<size_t>(a.size()); i++){
//        fmt::print("{}: {} =? {}\n",i, a.data()[i] , b.data()[i]);
        if(std::abs(a.data()[i] - b.data()[i])  > prec) return false;
    }
    return true;
}


TEST_CASE("Compare matrix and tensor kronecker products", "[kronecker]") {
    using cplx = std::complex<double>;
    SECTION("Rank 4") {
        Eigen::MatrixXcd id = qm::spinHalf::id;
        Eigen::MatrixXcd m(2,2);
        for(long i = 0; i < m.size(); i ++) m(i) = static_cast<double>(i+1);
        m(1,0) = -2;
        m(0,1) = std::complex<double>(-2.3,-10);

        auto matrix_m_id = linalg::matrix::kronecker(m,id);
        auto matrix_id_m = linalg::matrix::kronecker(id,m);
        auto tensor4_m_id = linalg::tensor::kronecker(m,id);
        auto tensor4_id_m = linalg::tensor::kronecker(id,m);
        fmt::print("matrix_m_id\n{}\n",  to_string(matrix_m_id));
        fmt::print("matrix_id_m\n{}\n",  to_string(matrix_id_m));
        fmt::print("tensor4_m_id\n{}\n", linalg::tensor::to_string(tensor4_m_id));
        fmt::print("tensor4_id_m\n{}\n", linalg::tensor::to_string(tensor4_id_m));
        REQUIRE(equal(matrix_m_id, tensor4_id_m));
        REQUIRE(equal(matrix_id_m, tensor4_m_id));
    }
    SECTION("Rank 6") {
        Eigen::MatrixXcd id = qm::spinHalf::id;
        Eigen::MatrixXcd m(2,2);
        for(long i = 0; i < m.size(); i ++) m(i) = static_cast<double>(i+1);
        m(1,0) = -2;
        m(0,1) = std::complex<double>(-2.3,-10);

        auto matrix_m_id_id = linalg::matrix::kronecker(linalg::matrix::kronecker(m,id),id);
        auto matrix_id_id_m = linalg::matrix::kronecker(linalg::matrix::kronecker(id,id),m);
        auto tensor6_m_id_id = linalg::tensor::kronecker(linalg::tensor::kronecker(m,id),id);
        auto tensor6_id_id_m = linalg::tensor::kronecker(linalg::tensor::kronecker(id,id),m);
        fmt::print("matrix_m_id_id\n{}\n", to_string(matrix_m_id_id));
        fmt::print("matrix_id_id_m\n{}\n", to_string(matrix_id_id_m));
        fmt::print("tensor6_m_id_id\n{}\n", linalg::tensor::to_string(tensor6_m_id_id));
        fmt::print("tensor6_id_id_m\n{}\n", linalg::tensor::to_string(tensor6_id_id_m));
        REQUIRE(equal(matrix_m_id_id, tensor6_id_id_m));
        REQUIRE(equal(matrix_id_id_m, tensor6_m_id_id));
    }
}


TEST_CASE("Test partial trace of tensors", "[partial trace]") {
    using cplx = std::complex<double>;
    SECTION("Rank 4") {
        Eigen::MatrixXcd id = qm::spinHalf::id;
        Eigen::MatrixXcd m(2,2);
        for(long i = 0; i < m.size(); i ++) m(i) = static_cast<double>(i+1);
        m(1,0) = -2;
        m(0,1) = std::complex<double>(-2.3,-10);

        // Compute the traces
        auto tr_id = id.trace();
        auto tr_m  = m.trace();
        // Define the expected results
        Eigen::MatrixXcd matrix_m = m * tr_id;
        Eigen::MatrixXcd matrix_id = id * tr_m;
        // Define the tensors to trace
        auto tensor4_m_id = linalg::tensor::kronecker(m,id);
        auto tensor4_id_m = linalg::tensor::kronecker(id,m);
        // Do the tracing
        auto tensor2_m = linalg::tensor::trace(tensor4_m_id,Textra::idx({1},{3})) ;
        auto tensor2_id = linalg::tensor::trace(tensor4_m_id,Textra::idx({0},{2}));
        fmt::print("matrix_m\n{}\n", to_string(matrix_m));
        fmt::print("matrix_id\n{}\n", to_string(matrix_id));
        fmt::print("tensor2_m\n{}\n", to_string(tensor2_m));
        fmt::print("tensor2_id\n{}\n", to_string(tensor2_id));

        REQUIRE(equal(matrix_m, tensor2_m));
        REQUIRE(equal(matrix_id, tensor2_id));

        // Now do the trace once more to get a scalar
        auto trace_m_a = linalg::tensor::trace(tensor2_m,Textra::idx({0},{1}));
        auto trace_id_a = linalg::tensor::trace(tensor2_id,Textra::idx({0},{1}));

        REQUIRE(trace_m_a(0) == tr_id * tr_m);
        REQUIRE(trace_id_a(0) ==  tr_id * tr_m);

        // Trace directly down to scalar
        auto trace_m_b = linalg::tensor::trace(tensor4_m_id,Textra::idx({0,1},{2,3}));
        auto trace_id_b = linalg::tensor::trace(tensor4_m_id,Textra::idx({0,1},{2,3}));

        REQUIRE(trace_m_b(0) == tr_id * tr_m);
        REQUIRE(trace_id_b(0) ==  tr_id * tr_m);
    }
    SECTION("Rank 6") {
        Eigen::MatrixXcd id = qm::spinHalf::id;
        Eigen::MatrixXcd m(2,2);
        for(long i = 0; i < m.size(); i ++) m(i) = static_cast<double>(i+1);
        m(1,0) = -2;
        m(0,1) = std::complex<double>(-2.3,-10);

        // Compute the traces
        auto tr_id = id.trace();
        auto tr_m  = m.trace();
        // Define the expected results
        Eigen::MatrixXcd matrix_m = m * tr_id * tr_id;
        Eigen::MatrixXcd matrix_id = id * tr_m * tr_id;
        Eigen::MatrixXcd matrix_m_id = linalg::matrix::kronecker(id,m) * tr_id;
        Eigen::MatrixXcd matrix_id_m = linalg::matrix::kronecker(m,id) * tr_id;

        // Define the tensors to trace
        auto tensor6_m_id_id = linalg::tensor::kronecker(linalg::tensor::kronecker(m,id),id);
        auto tensor6_id_m_id = linalg::tensor::kronecker(linalg::tensor::kronecker(id,m),id);
        auto tensor6_id_id_m = linalg::tensor::kronecker(linalg::tensor::kronecker(id,id),m);
        // Do the tracing one step
        auto tensor4_m_id_a = linalg::tensor::trace(tensor6_m_id_id,Textra::idx({2},{5}));
        auto tensor4_m_id_b = linalg::tensor::trace(tensor6_id_m_id,Textra::idx({0},{3}));
        auto tensor4_id_m_a = linalg::tensor::trace(tensor6_id_m_id,Textra::idx({2},{5}));
        auto tensor4_id_m_b = linalg::tensor::trace(tensor6_id_id_m,Textra::idx({0},{3}));

        REQUIRE(equal(tensor4_m_id_a, tensor4_m_id_b));
        REQUIRE(equal(tensor4_id_m_a, tensor4_id_m_b));
        REQUIRE(equal(tensor4_m_id_a, matrix_m_id));
        REQUIRE(equal(tensor4_id_m_a, matrix_id_m));

        // Do the tracing once more
        auto tensor2_m_a = linalg::tensor::trace(tensor4_m_id_a,Textra::idx({1},{3}));
        auto tensor2_m_b = linalg::tensor::trace(tensor4_id_m_a,Textra::idx({0},{2}));
        auto tensor2_id_a = linalg::tensor::trace(tensor4_m_id_a,Textra::idx({0},{2}));
        auto tensor2_id_b = linalg::tensor::trace(tensor4_id_m_a,Textra::idx({1},{3}));

        REQUIRE(equal(tensor2_m_a, tensor2_m_b));
        REQUIRE(equal(tensor2_id_a, tensor2_id_a));
        REQUIRE(equal(tensor2_m_a, matrix_m));
        REQUIRE(equal(tensor2_id_a, matrix_id));

        // Do the tracing in two steps directly
        auto tensor2_m_c = linalg::tensor::trace(tensor6_m_id_id,Textra::idx({1,2},{4,5}));
        auto tensor2_m_d = linalg::tensor::trace(tensor6_id_m_id,Textra::idx({0,2},{3,5}));
        auto tensor2_m_e = linalg::tensor::trace(tensor6_id_id_m,Textra::idx({0,1},{3,4}));

        REQUIRE(equal(tensor2_m_a, tensor2_m_c));
        REQUIRE(equal(tensor2_m_a, tensor2_m_d));
        REQUIRE(equal(tensor2_m_a, tensor2_m_e));
    }
}



int main(int argc, char **argv){
    return Catch::Session().run(argc, argv);
}