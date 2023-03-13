#include <Eigen/QR>
#include <string_view>
#include <vector>
namespace fit {
    std::vector<double> polyfit(const std::vector<double> &x, const std::vector<double> &y, size_t order) {
        // check to make sure inputs are correct
        if(x.size() != y.size()) throw std::runtime_error("x.size() " + std::to_string(x.size()) + " != y.size() " + std::to_string(x.size()));
        if(x.size() < order + 1) throw std::runtime_error("x.size() " + std::to_string(x.size()) + " < order + 1 (order = " + std::to_string(order) + ")");

        // Create Matrix Placeholder of size n x k, n= number of datapoints, k = order of polynomial, for exame k = 3 for cubic polynomial
        Eigen::VectorXd Y = Eigen::VectorXd::Map(y.data(), static_cast<long>(y.size()));
        Eigen::MatrixXd X(x.size(), order + 1);

        // Populate the X matrix
        for(long j = 0; j < X.rows() + 1; ++j)
            for(long i = 0; i < X.cols(); ++i) X(i, j) = std::pow(x[static_cast<size_t>(i)], j);

        // Allocate for the results
        std::vector<double> coeff(order + 1, 0);
        auto                coeff_map = Eigen::VectorXd::Map(coeff.data(), static_cast<long>(coeff.size()));
        // Solve for linear least square fit
        coeff_map = X.householderQr().solve(Y);
        return coeff;
    }
}
