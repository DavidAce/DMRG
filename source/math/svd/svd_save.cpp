//
// Created by david on 2021-11-12.
//

#include "../svd.h"
#include <h5pp/h5pp.h>

template<typename Scalar>
void svd::solver::save_svd(const MatrixType<Scalar> &A, const MatrixType<Scalar> &U, const VectorType<Scalar> &S, const MatrixType<Scalar> &VT, long rank_max,
                           const std::string &lib, const std::vector<std::pair<std::string, std::string>> &details) const {
    auto rank = S.size();
    auto rows = A.rows();
    auto cols = A.cols();

    auto file       = h5pp::File("svd-save.h5", h5pp::FilePermission::READWRITE);
    auto group_num  = 0;
    auto group_name = fmt::format("svd_{}_{}", lib, group_num);
    while(file.linkExists(group_name)) group_name = fmt::format("svd_{}_{}", lib, ++group_num);
    file.writeDataset(A, fmt::format("{}/A", group_name));
    file.writeDataset(U, fmt::format("{}/U", group_name));
    file.writeDataset(S, fmt::format("{}/S", group_name));
    file.writeDataset(VT, fmt::format("{}/VT", group_name));
    file.writeAttribute(rows, "rows", group_name);
    file.writeAttribute(cols, "cols", group_name);
    file.writeAttribute(rank, "rank", group_name);
    file.writeAttribute(rank_max, "rank_max", group_name);
    file.writeAttribute(use_bdc, "use_bdc", group_name);
    file.writeAttribute(threshold, "threshold", group_name);
    file.writeAttribute(switchsize_bdc, "switchsize_bdc", group_name);
    for(const auto &[key, val] : details) file.writeAttribute(val, key, group_name);

    //
    //
    //
    //#if defined(OPENBLAS_AVAILABLE)
    //    file.writeAttribute(OPENBLAS_VERSION, "OPENBLAS_VERSION", group_name);
    //    file.writeAttribute(openblas_get_num_threads(), "openblas_get_num_threads", group_name);
    //    file.writeAttribute(openblas_get_parallel(), "openblas_parallel_mode", group_name);
    //    file.writeAttribute(openblas_get_corename(), "openblas_get_corename", group_name);
    //    file.writeAttribute(openblas_get_config(), "openblas_get_config()", group_name);
    //    file.writeAttribute(OPENBLAS_GEMM_MULTITHREAD_THRESHOLD, "OPENBLAS_GEMM_MULTITHREAD_THRESHOLD", group_name);
    //#endif
    //
    //#if defined(MKL_AVAILABLE)
    //    MKLVersion Version;
    //    mkl_get_version(&Version);
    //    file.writeAttribute(Version.MajorVersion, "Intel-MKL-MajorVersion", group_name);
    //    file.writeAttribute(Version.MinorVersion, "Intel-MKL-MinorVersion", group_name);
    //    file.writeAttribute(Version.UpdateVersion, "Intel-MKL-UpdateVersion", group_name);
    //#endif
}
using cplx = std::complex<double>;
using real = double;

template void svd::solver::save_svd(const MatrixType<real> &A, const MatrixType<real> &U, const VectorType<real> &S, const MatrixType<real> &VT, long rank_max,
                                    const std::string &lib, const std::vector<std::pair<std::string, std::string>> &details) const;
template void svd::solver::save_svd(const MatrixType<cplx> &A, const MatrixType<cplx> &U, const VectorType<cplx> &S, const MatrixType<cplx> &VT, long rank_max,
                                    const std::string &lib, const std::vector<std::pair<std::string, std::string>> &details) const;