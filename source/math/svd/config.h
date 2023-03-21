#pragma once
#include <fmt/core.h>
#include <fmt/std.h>
#include <optional>
#include <string>
#include <variant>

namespace svd {
    enum class lib { eigen, lapacke }; /*!< Library */

    /*! \brief SVD routine. In Eigen gesvd, gejsv, gesvj, are just JacobiSVD, while gesdd is BDSVD
     * Links to disussions:
     * https://discourse.julialang.org/t/svd-better-default-to-gesvd-instead-of-gesdd/20603/17
     * https://mathoverflow.net/questions/161252/what-is-the-time-complexity-of-truncated-svd
     *
     *
     * */
    enum class rtn {
        gesvj,  /*!< Slowest. Preconditioned Jacobi. Probably the most accurate for tiny singular values  */
        gejsv,  /*!< Slower. Preconditioned Jacobi. More accurate than QR */
        gesvd,  /*!< Fast. Bidiagonal QR iteration. Sufficient accuracy. */
        gesvdx, /*!< Faster for few k, otherwise like gesvd. Computes the top k singular values by calling an eigenvalue solver on [0,A; A^T,0]   */
        gesdd,  /*!< Fastest for all k. Scalable Divide and conquer approach with Bidiagonal QR preconditioner?. Low accuracy on small singular values.   */
        gersvd, /*!< Fastest for very few k. Randomized SVD works well for very low rank approximations of huge matrices */
        geauto  /*!< Defaults to gejsv for small matrices, moving to gesdd for large matrices, and to either gesvdx or gersvd for low k */
    };

    using svdx_select_t = std::variant<size_t, double>;

    constexpr inline std::string_view enum2sv(svd::lib lib) {
        switch(lib) {
            case svd::lib::eigen: return "eigen";
            case svd::lib::lapacke: return "lapacke";
            default: throw std::logic_error("Could not match svd::lib");
        }
    }

    constexpr inline std::string_view enum2sv(svd::rtn rtn) {
        switch(rtn) {
            case svd::rtn::gesvd: return "gesvd";
            case svd::rtn::gejsv: return "gejsv";
            case svd::rtn::gesvj: return "gesvj";
            case svd::rtn::gesdd: return "gesdd";
            case svd::rtn::gesvdx: return "gesvdx";
            case svd::rtn::gersvd: return "gersvd";
            case svd::rtn::geauto: return "geauto";
            default: throw std::logic_error("Could not match svd::rtn");
        }
    }

    struct config {
        std::optional<long>          rank_max         = std::nullopt;
        std::optional<long>          rank_min         = std::nullopt; /*!< Keep this many singular values even if they are smaller than truncation_lim */
        std::optional<double>        truncation_limit = std::nullopt;
        std::optional<size_t>        switchsize_gejsv = std::nullopt;
        std::optional<size_t>        switchsize_gesvd = std::nullopt;
        std::optional<size_t>        switchsize_gesdd = std::nullopt;
        std::optional<svdx_select_t> svdx_select      = std::nullopt;
        std::optional<size_t>        loglevel         = std::nullopt;
        std::optional<svd::lib>      svd_lib          = std::nullopt;
        std::optional<svd::rtn>      svd_rtn          = std::nullopt;
        std::optional<bool>          save_fail        = std::nullopt;
        std::optional<bool>          benchmark        = std::nullopt;
        std::string                  to_string();
        config() = default;
        explicit config(long rank_max_);
        explicit config(double truncation_lim_);
        explicit config(long rank_max_, double truncation_lim_);
        explicit config(std::optional<long> rank_max_);
        explicit config(std::optional<double> truncation_lim_);
        explicit config(std::optional<long> rank_max_, std::optional<double> truncation_lim_);
    };
}