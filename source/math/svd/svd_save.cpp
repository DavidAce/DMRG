//
// Created by david on 2021-11-12.
//

#if defined(MKL_AVAILABLE)
    #ifndef SVD_SAVE_OPENBLAS_ATTRIBUTES
        #define SVD_SAVE_MKL_ATTRIBUTES
    #endif
    #include <mkl.h>
#elif defined(OPENBLAS_AVAILABLE)
    #ifndef SVD_SAVE_OPENBLAS_ATTRIBUTES
        #define SVD_SAVE_OPENBLAS_ATTRIBUTES
    #endif
    #include <openblas/cblas.h>
#elif defined(FLEXIBLAS_AVAILABLE)
    #include <flexiblas/flexiblas_config.h>
    #ifndef SVD_SAVE_FLEXIBLAS_ATTRIBUTES
        #define SVD_SAVE_FLEXIBLAS_ATTRIBUTES
    #endif
#elif __has_include(<cblas-openblas.h>)
    #ifndef SVD_SAVE_OPENBLAS_ATTRIBUTES
        #define SVD_SAVE_OPENBLAS_ATTRIBUTES
    #endif
    #include <cblas-openblas.h>
#endif

#include "../svd.h"
#include "config/settings.h"
#include <Eigen/src/Core/util/Macros.h>
#include <h5pp/h5pp.h>

void svd::solver::save_svd() {
    auto &smd = saveMetaData;
    if(smd.svd_save == save::NONE) return;
    if(not smd.svd_is_running) return;
    auto directory = h5pp::fs::path(settings::storage::output_filepath).parent_path().string();
    auto filepath  = fmt::format("{}/svd-save-{}.h5", directory, settings::input::seed);
    svd::log->info("Saving SVD to file: {}", filepath);
    auto file       = h5pp::File(filepath, h5pp::FilePermission::READWRITE);
    auto group_num  = 0;
    auto group_name = fmt::format("svd_{}", group_num);
    if(smd.svd_save == save::ALL)
        while(file.linkExists(group_name)) group_name = fmt::format("svd_{}", ++group_num);
    if(smd.svd_save == save::LAST) group_name = "svd-last";
    if(smd.svd_save == save::FAIL) group_name = "svd-fail";
    using MatrixReal = svd::internal::SaveMetaData::MatrixReal;
    using MatrixCplx = svd::internal::SaveMetaData::MatrixCplx;
    if(auto *A = std::get_if<MatrixReal>(&smd.A); A != nullptr and A->size() != 0) {
        file.writeDataset(*A, fmt::format("{}/A_real", group_name), H5D_layout_t::H5D_CHUNKED);
    }
    if(auto *A = std::get_if<MatrixCplx>(&smd.A); A != nullptr and A->size() != 0) {
        file.writeDataset(*A, fmt::format("{}/A_cplx", group_name), H5D_layout_t::H5D_CHUNKED);
    }
    if(auto *U = std::get_if<MatrixReal>(&smd.U); U != nullptr and U->size() != 0) {
        file.writeDataset(*U, fmt::format("{}/U_real", group_name), H5D_layout_t::H5D_CHUNKED);
    }
    if(auto *U = std::get_if<MatrixCplx>(&smd.U); U != nullptr and U->size() != 0) {
        file.writeDataset(*U, fmt::format("{}/U_cplx", group_name), H5D_layout_t::H5D_CHUNKED);
    }
    if(smd.S.size() != 0) { file.writeDataset(smd.S, fmt::format("{}/S", group_name), H5D_layout_t::H5D_CHUNKED); }
    if(auto *VT = std::get_if<MatrixReal>(&smd.VT); VT != nullptr and VT->size() != 0) {
        file.writeDataset(*VT, fmt::format("{}/VT_real", group_name), H5D_layout_t::H5D_CHUNKED);
    }
    if(auto *VT = std::get_if<MatrixCplx>(&smd.VT); VT != nullptr and VT->size() != 0) {
        file.writeDataset(*VT, fmt::format("{}/VT_cplx", group_name), H5D_layout_t::H5D_CHUNKED);
    }

    file.writeAttribute(settings::input::seed, group_name, "seed");
    file.writeAttribute(smd.rank_max, group_name, "rank_max");
    file.writeAttribute(enum2sv(smd.svd_lib), group_name, "svd_lib");
    file.writeAttribute(enum2sv(smd.svd_rtn), group_name, "svd_rtn");
    file.writeAttribute(smd.truncation_lim, group_name, "truncation_lim");
    file.writeAttribute(smd.switchsize_gejsv, group_name, "switchsize_gejsv");
    file.writeAttribute(smd.switchsize_gesvd, group_name, "switchsize_gesvd");
    file.writeAttribute(smd.switchsize_gesdd, group_name, "switchsize_gesdd");
    file.writeAttribute(smd.info, group_name, "info");
    file.writeAttribute(smd.truncation_error, group_name, "truncation_error");

    if(smd.svd_lib == svd::lib::lapacke) {
#if defined(SVD_SAVE_OPENBLAS_ATTRIBUTES)
        file.writeAttribute(OPENBLAS_VERSION, "OPENBLAS_VERSION", group_name);
        file.writeAttribute(openblas_get_num_threads(), "openblas_get_num_threads", group_name);
        file.writeAttribute(openblas_get_parallel(), "openblas_parallel_mode", group_name);
        file.writeAttribute(openblas_get_corename(), "openblas_get_corename", group_name);
        file.writeAttribute(openblas_get_config(), "openblas_get_config()", group_name);
        file.writeAttribute(OPENBLAS_GEMM_MULTITHREAD_THRESHOLD, "OPENBLAS_GEMM_MULTITHREAD_THRESHOLD", group_name);
#endif

#if defined(SVD_SAVE_MKL_ATTRIBUTES)
        MKLVersion Version;
        mkl_get_version(&Version);
        file.writeAttribute(Version.MajorVersion, "Intel-MKL-MajorVersion", group_name);
        file.writeAttribute(Version.MinorVersion, "Intel-MKL-MinorVersion", group_name);
        file.writeAttribute(Version.UpdateVersion, "Intel-MKL-UpdateVersion", group_name);
#endif
#if defined(SVD_SAVE_FLEXIBLAS_ATTRIBUTES)
        file.writeAttribute(FLEXIBLAS_DEFAULT_LIB_PATH, group_name, "FLEXIBLAS_DEFAULT_LIB_PATH");
        file.writeAttribute(FLEXIBLAS_VERSION, group_name, "FLEXIBLAS_VERSION");
#endif

    } else if(smd.svd_lib == svd::lib::eigen) {
        auto eigen_version = fmt::format("{}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
        file.writeAttribute(eigen_version, group_name, "Eigen Version");
    }
    smd = svd::internal::SaveMetaData{};
//    if (group_num > 50) std::exit(0);
}

template<typename Scalar>
void svd::solver::save_svd(const MatrixType<Scalar> &A, std::optional<svd::save> override) const {
    svd::save save_internal = override.has_value() ? override.value() : svd_save;
    if(save_internal == save::NONE) return;
    if(save_internal == save::FAIL) return;

    auto rows      = A.rows();
    auto cols      = A.cols();
    auto directory = h5pp::fs::path(settings::storage::output_filepath).parent_path().string();
    auto filepath  = fmt::format("{}/svd-save-{}.h5", directory, settings::input::seed);
    svd::log->debug("Saving A to file: {}", filepath);
    auto file       = h5pp::File(filepath, h5pp::FilePermission::READWRITE);
    auto group_num  = 0;
    auto group_name = fmt::format("svd_{}", group_num);
    if(override) {}
    if(save_internal == save::ALL or override == save::ALL)
        while(file.linkExists(group_name)) group_name = fmt::format("svd_{}", ++group_num);
    if(save_internal == save::LAST or override == save::LAST) group_name = "svd-last";
    file.writeDataset(A, fmt::format("{}/A", group_name), H5D_layout_t::H5D_CHUNKED);
    file.writeAttribute(rows, group_name, "rows");
    file.writeAttribute(cols, group_name, "cols");
    file.writeAttribute(settings::input::seed, group_name, "seed");
    file.writeAttribute(rank_max, group_name, "rank_max");
    file.writeAttribute(enum2sv(svd_lib), group_name, "svd_lib");
    file.writeAttribute(enum2sv(svd_rtn), group_name, "svd_rtn");
    file.writeAttribute(truncation_lim, group_name, "truncation_lim");
    file.writeAttribute(switchsize_gejsv, group_name, "switchsize_gejsv");
    file.writeAttribute(switchsize_gesvd, group_name, "switchsize_gesvd");
    file.writeAttribute(switchsize_gesdd, group_name, "switchsize_gesdd");

    if(svd_lib == svd::lib::lapacke) {
#if defined(SVD_SAVE_OPENBLAS_ATTRIBUTES)
        file.writeAttribute(OPENBLAS_VERSION, "OPENBLAS_VERSION", group_name);
        file.writeAttribute(openblas_get_num_threads(), "openblas_get_num_threads", group_name);
        file.writeAttribute(openblas_get_parallel(), "openblas_parallel_mode", group_name);
        file.writeAttribute(openblas_get_corename(), "openblas_get_corename", group_name);
        file.writeAttribute(openblas_get_config(), "openblas_get_config()", group_name);
        file.writeAttribute(OPENBLAS_GEMM_MULTITHREAD_THRESHOLD, "OPENBLAS_GEMM_MULTITHREAD_THRESHOLD", group_name);
#endif

#if defined(SVD_SAVE_MKL_ATTRIBUTES)
        MKLVersion Version;
        mkl_get_version(&Version);
        file.writeAttribute(Version.MajorVersion, "Intel-MKL-MajorVersion", group_name);
        file.writeAttribute(Version.MinorVersion, "Intel-MKL-MinorVersion", group_name);
        file.writeAttribute(Version.UpdateVersion, "Intel-MKL-UpdateVersion", group_name);
#endif
#if defined(SVD_SAVE_FLEXIBLAS_ATTRIBUTES)
        file.writeAttribute(FLEXIBLAS_DEFAULT_LIB_PATH, group_name, "FLEXIBLAS_DEFAULT_LIB_PATH");
        file.writeAttribute(FLEXIBLAS_VERSION, group_name, "FLEXIBLAS_VERSION");
#endif
    } else if(svd_lib == svd::lib::eigen) {
        auto eigen_version = fmt::format("{}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
        file.writeAttribute(eigen_version, group_name, "Eigen Version");
    }
}

template<typename Scalar>
void svd::solver::save_svd(const MatrixType<Scalar> &U, const VectorType<Scalar> &S, const MatrixType<Scalar> &VT, int info) const {
    if(svd_save == save::NONE) return;
    if(svd_save == save::FAIL) return;
    auto directory  = h5pp::fs::path(settings::storage::output_filepath).parent_path().string();
    auto filepath   = fmt::format("{}/svd-save-{}.h5", directory, settings::input::seed);
    auto file       = h5pp::File(filepath, h5pp::FilePermission::READWRITE);
    auto group_num  = 0;
    auto group_name = fmt::format("svd_{}", group_num);

    if(svd_save == save::ALL)
        while(file.linkExists(group_name)) group_name = fmt::format("svd_{}", ++group_num);
    if(svd_save == save::LAST) group_name = "svd-last";
    file.writeDataset(U, fmt::format("{}/U", group_name), H5D_layout_t::H5D_CHUNKED);
    file.writeDataset(S, fmt::format("{}/S", group_name), H5D_layout_t::H5D_CHUNKED);
    file.writeDataset(VT, fmt::format("{}/VT", group_name), H5D_layout_t::H5D_CHUNKED);
    file.writeAttribute(S.size(), group_name, "rank");
    file.writeAttribute(info, group_name, "info");
}

using cplx = std::complex<double>;
using real = double;
template void svd::solver::save_svd(const MatrixType<real> &A, std::optional<svd::save> override) const;
template void svd::solver::save_svd(const MatrixType<cplx> &A, std::optional<svd::save> override) const;
template void svd::solver::save_svd(const MatrixType<real> &U, const VectorType<real> &S, const MatrixType<real> &VT, int info) const;
template void svd::solver::save_svd(const MatrixType<cplx> &U, const VectorType<cplx> &S, const MatrixType<cplx> &VT, int info) const;