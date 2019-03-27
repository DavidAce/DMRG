//
// Created by david on 2018-10-29.
//
//#ifndef HAVE_LAPACK_CONFIG_H
//#define HAVE_LAPACK_CONFIG_H
//#endif


#include <complex>

#ifdef MKL_AVAILABLE
#ifndef MKL_Complex8
#define MKL_Complex8 std::complex<float>
#endif
#ifndef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#endif
#ifndef ComplexFloat_
#define ComplexFloat_  MKL_Complex8
#endif
#ifndef ComplexDouble_
#define ComplexDouble_ MKL_Complex16
#endif
#include <mkl_lapacke.h>
#else
#define ComplexFloat_  __complex__ float
#define ComplexDouble_ __complex__ double
#include <lapacke.h>
#endif
#include "class_eigsolver.h"





class_eigsolver::class_eigsolver(size_t logLevel_):logLevel(logLevel_){
    eigutils::eigLogger::setLogger("eig",logLevel,false);
}


void class_eigsolver::setLogLevel(size_t logLevelZeroToSix){
    logLevel = logLevelZeroToSix;
    eigutils::eigLogger::setLogLevel(logLevelZeroToSix);
}


void class_eigsolver::eigs_init(const int L,
                                const int nev,
                                const int ncv,
                                const std::complex<double> sigma,
                                const eigutils::eigSetting::Type type,
                                const eigutils::eigSetting::Form form,
                                const eigutils::eigSetting::Ritz ritz,
                                const eigutils::eigSetting::Side side,
                                const eigutils::eigSetting::Storage storage,
                                const bool compute_eigvecs_,
                                const bool remove_phase_
)
{
    eigutils::eigLogger::log->trace("Init full");
    using namespace eigutils::eigSetting;
    solution.reset();
    bool is_shifted             = sigma == sigma;
    solverConf.compute_eigvecs  = compute_eigvecs_;
    solverConf.remove_phase     = remove_phase_;
    solverConf.eigMaxNev        = nev;
    solverConf.eigMaxNcv        = ncv;
    solverConf.shift            = is_shifted ? Shift::ON : Shift::OFF;
    solverConf.sigma            = sigma;
    solverConf.type             = type;
    solverConf.form             = form;
    solverConf.ritz             = ritz;
    solverConf.side             = side;
    solverConf.storage          = storage;


    if (ncv < nev ){
        if (nev >= 1 and nev <= 16 ){
            solverConf.eigMaxNcv = 8 + std::ceil((int)(1.5*nev));
        }
        else
        if (nev > 16 and nev <= L ){
            solverConf.eigMaxNcv = 2*nev;
        }
    }

    if (solverConf.form == Form::NONSYMMETRIC){
        if (solverConf.eigMaxNev == 1) {
            solverConf.eigMaxNev = 2;
        }
    }
    if (solverConf.eigMaxNcv >= L ){
        solverConf.eigMaxNcv = (L + nev)  / 2;
    }

    assert(solverConf.eigMaxNcv <= L and "Ncv > L");
    assert(solverConf.eigMaxNcv >= solverConf.eigMaxNev and "Ncv < Nev");
    assert(solverConf.eigMaxNev <= L and "Nev > L");
    solverConf.confOK = true;
}


void class_eigsolver::eig_init  ( const eigutils::eigSetting::Type type      ,
                                  const eigutils::eigSetting::Form form      ,
                                  const eigutils::eigSetting::Side side      ,
                                  const bool compute_eigvecs_                ,
                                  const bool remove_phase_
)
{
    eigutils::eigLogger::log->trace("Init light");
    using namespace eigutils::eigSetting;
    solution.reset();
    bool is_shifted             = false;
    solverConf.compute_eigvecs  = compute_eigvecs_;
    solverConf.remove_phase     = remove_phase_;
    solverConf.eigMaxNev        = 0;
    solverConf.eigMaxNcv        = 0;
    solverConf.shift            = is_shifted ? Shift::ON : Shift::OFF;
    solverConf.sigma            = std::numeric_limits<double>::quiet_NaN();
    solverConf.type             = type;
    solverConf.form             = form;
    solverConf.ritz             = Ritz::LM;
    solverConf.side             = side;
    solverConf.storage          = Storage::DENSE;
    solverConf.confOK           = true;
}




int class_eigsolver::eig_dsyevd(const double* matrix, int L){
    eigutils::eigLogger::log->trace("Starting eig_dsyevd. Eigvecs: {}", solverConf.compute_eigvecs);
    using namespace eigutils::eigSetting;
    auto & eigvals = solution.get_eigvals<Form::SYMMETRIC>();
    auto & eigvecs = solution.get_eigvecs<Type::REAL,Form::SYMMETRIC>();
    eigvals.resize(L);
    eigvecs.resize(L*L);
    std::copy(matrix, matrix + L*L, eigvecs.begin());
    int info = eig_dsyevd(eigvecs.data(),eigvals.data(), L);
    if (info == 0){
        solution.meta.eigvecs_found = true;
        solution.meta.eigvals_found = true;
        solution.meta.rows           = L;
        solution.meta.cols           = L;
        solution.meta.nev            = L;
        solution.meta.n              = L;
        solution.meta.form           = Form::SYMMETRIC;
        solution.meta.type           = Type::REAL ;
    }else{
        throw std::runtime_error("LAPACK dsyevd failed with error: " + std::to_string(info));
    }
    return info;
}

int class_eigsolver::eig_dsyevd(double *matrix2eigvecs, double * eigvals, int L){
    eigutils::eigLogger::log->trace("Starting eig_dsyevd (lapacke). Eigvecs: {}", solverConf.compute_eigvecs);


    // For some reason the recommended lwork from netlib doesn't work. It's better to ask lapack with a query.
    // These nice values are inspired from armadillo. The prefactors give good performance.
//    int lwork  = 2 * (1 + 6*L + 2*(L*L));
//    int liwork = 3 * (3 + 5*L);
    int info   = 0;
    char jobz = solverConf.compute_eigvecs ? 'V' : 'N';
    double lwork_query [1];
    int    liwork_query[1];

    info = LAPACKE_dsyevd_work(LAPACK_COL_MAJOR,jobz,'U',L,
                               matrix2eigvecs,
                               L,
                               eigvals,
                               lwork_query,
                               -1,
                               liwork_query,
                               -1);

    int lwork     = (int) 2 * lwork_query[0]; //Make it twice as big for performance.
    int liwork    = (int) 3 * liwork_query[0]; //Make it thrice as big for performance.
    eigutils::eigLogger::log->trace(" lwork  = {}", lwork);
    eigutils::eigLogger::log->trace(" liwork = {}", liwork);

    std::vector<double> work  ( (unsigned long) lwork );
    std::vector<int   > iwork ( (unsigned long) liwork );

    info = LAPACKE_dsyevd_work(LAPACK_COL_MAJOR,jobz,'U',L,
                               matrix2eigvecs,
                               L,
                               eigvals,
                               work.data(),
                               lwork,
                               iwork.data(),
                               liwork);
    return info;
}



int class_eigsolver::eig_zheevd(const std::complex<double>* matrix, int L){
    eigutils::eigLogger::log->trace("Starting eig_zheevd");
    using namespace eigutils::eigSetting;
    auto & eigvals = solution.get_eigvals<Form::SYMMETRIC>();
    auto & eigvecs = solution.get_eigvecs<Type::CPLX,Form::SYMMETRIC>();
    eigvals.resize(L);
    eigvecs.resize(L*L);
    std::copy(matrix, matrix + L*L, eigvecs.begin());
    int info = eig_zheevd(eigvecs.data(), eigvals.data(), L);
    if (info == 0){
        solution.meta.eigvecs_found = true;
        solution.meta.eigvals_found = true;
        solution.meta.rows           = L;
        solution.meta.cols           = L;
        solution.meta.nev            = L;
        solution.meta.n              = L;
        solution.meta.form           = Form::SYMMETRIC;
        solution.meta.type           = Type::CPLX ;
    }else{
        throw std::runtime_error("LAPACK zheevd failed with error: " + std::to_string(info));
    }
    return info;
}

int class_eigsolver::eig_zheevd(std::complex<double>* matrix2eigvecs, double *eigvals, int L){
    eigutils::eigLogger::log->trace("Starting eig_zheevd (lapacke)");
    using Scalar = std::complex<double>;
    //These nice values are inspired from armadillo. The prefactors give good performance.
    int lwork  = 2 * (2*L + L*L);
    int lrwork = 2 * (1 + 5*L + 2*(L*L));
    int liwork = 3 * (3 + 5*L);
    int info   = 0;
    std::vector<Scalar> work  ( lwork );
    std::vector<double> rwork ( lrwork );
    std::vector<int   > iwork ( liwork );
    char jobz = solverConf.compute_eigvecs ? 'V' : 'N';
    info = LAPACKE_zheevd_work(LAPACK_COL_MAJOR,jobz,'U',L,
            reinterpret_cast< ComplexDouble_ *>(matrix2eigvecs),
            L,
            eigvals,
            reinterpret_cast< ComplexDouble_ *>(work.data()),
            lwork,
            rwork.data(),
            lrwork,
            iwork.data(),
            liwork);
    return info;
}


int class_eigsolver::eig_dgeev(const double* matrix, int L){
    eigutils::eigLogger::log->trace("Starting eig_dgeev");
    using namespace eigutils::eigSetting;
    auto & eigvals  = solution.get_eigvals<Form::NONSYMMETRIC>();
    auto & eigvecsR = solution.get_eigvecs<Type::REAL,Form::NONSYMMETRIC, Side::R>();
    auto & eigvecsL = solution.get_eigvecs<Type::REAL,Form::NONSYMMETRIC, Side::L>();
    eigvals.resize(L);
    eigvecsR.resize(L*L);
    eigvecsL.resize(L*L);

    int info = eig_dgeev(matrix,eigvecsR.data(),eigvecsL.data(),eigvals.data(),L );
    if (info == 0){
        solution.meta.eigvecsR_found = true;
        solution.meta.eigvecsL_found = true;
        solution.meta.eigvals_found  = true;
        solution.meta.rows           = L;
        solution.meta.cols           = L;
        solution.meta.nev            = L;
        solution.meta.n              = L;

        solution.meta.form = Form::NONSYMMETRIC;
        solution.meta.type = Type::REAL ;
    }else{
        throw std::runtime_error("LAPACK dgeev failed with error: " + std::to_string(info));
    }
    return info;
}


int class_eigsolver::eig_dgeev(const double* matrix, std::complex<double> *eigvecsR, std::complex<double>* eigvecsL, std::complex<double> *eigvals, int L){
    eigutils::eigLogger::log->trace("Starting eig_dgeev (lapacke)");

    // For some reason the recommended lwork from netlib doesn't work. It's better to ask lapack with a query.
//    int lwork   = 2 * (4*L);
    int info    = 0;
    double lwork_query;

    std::vector<double> eigvals_real(L);
    std::vector<double> eigvals_imag(L);
    Eigen::MatrixXd tmpR(L,L);
    Eigen::MatrixXd tmpL(L,L);

    char jobz = solverConf.compute_eigvecs ? 'V' : 'N';
    info = LAPACKE_dgeev_work(LAPACK_COL_MAJOR,jobz,jobz,L,
                              const_cast<double*>(matrix),
                              L,
                              eigvals_real.data(),
                              eigvals_imag.data(),
                              tmpL.data(),
                              L,
                              tmpR.data(),
                              L,
                              &lwork_query,
                              -1);
    int lwork = (int) std::real(2.0*lwork_query); //Make it twice as big for performance.
    std::vector<double> work  ( (unsigned long)lwork );
    info = LAPACKE_dgeev_work(LAPACK_COL_MAJOR,jobz,jobz,L,
                              const_cast<double*>(matrix),
                              L,
                              eigvals_real.data(),
                              eigvals_imag.data(),
                              tmpL.data(),
                              L,
                              tmpR.data(),
                              L,
                              work.data(),
                              lwork);


    int count = 0;
    for (int i = 0; i < L; i++) {
        eigvals[i] = std::complex<double>(eigvals_real[i], eigvals_imag[i]);
        int j = 0;
        while (j < L){
            if (eigvals_imag[j] == 0.0){
                eigvecsR[count] = tmpR(i,j);//tmpR[i + j*L];
                eigvecsL[count] = tmpL(i,j);//tmpL[i + j*L];
                count++;
                j++;
            }else{
                eigvecsR[count] = std::complex<double>(tmpR(i,j),tmpR(i,j+1));   // std::complex<double>(tmpR[i + j*L], tmpR[i + (j+1)*L]);
                eigvecsL[count] = std::complex<double>(tmpL(i,j),tmpL(i,j+1));   //std::complex<double>(tmpL[i + j*L], tmpL[i + (j+1)*L]);
                count++;
                eigvecsR[count] = std::complex<double>(tmpR(i,j), -tmpR(i,j+1)); //std::complex<double>(tmpR[i + j*L], -tmpR[i + (j+1)*L]);
                eigvecsL[count] = std::complex<double>(tmpL(i,j), -tmpL(i,j+1)); //std::complex<double>(tmpL[i + j*L], -tmpL[i + (j+1)*L]);
                count++;
                j+=2;
            }
        }

    }
    return info;
}



int class_eigsolver::eig_zgeev(const std::complex<double>* matrix, int L){
    eigutils::eigLogger::log->trace("Starting eig_zgeev");
    using namespace eigutils::eigSetting;
    auto & eigvals  = solution.get_eigvals<Form::NONSYMMETRIC>();
    auto & eigvecsR = solution.get_eigvecs<Type::CPLX,Form::NONSYMMETRIC, Side::R>();
    auto & eigvecsL = solution.get_eigvecs<Type::CPLX,Form::NONSYMMETRIC, Side::L>();
    eigvals.resize(L);
    eigvecsR.resize(L*L);
    eigvecsL.resize(L*L);
    int info = eig_zgeev(matrix,eigvecsR.data(), eigvecsL.data(),eigvals.data(),L);
    if (info == 0){
        solution.meta.eigvecsR_found = true;
        solution.meta.eigvecsL_found = true;
        solution.meta.eigvals_found  = true;
        solution.meta.rows           = L;
        solution.meta.cols           = L;
        solution.meta.nev            = L;
        solution.meta.n              = L;
        solution.meta.form           = Form::NONSYMMETRIC;
        solution.meta.type           = Type::CPLX ;
    }else{
        throw std::runtime_error("LAPACK zgeev failed with error: " + std::to_string(info));
    }
    return info;
}




int class_eigsolver::eig_zgeev(const std::complex<double>* matrix, std::complex<double>* eigvecsR, std::complex<double>* eigvecsL, std::complex<double> *eigvals, int L){
    eigutils::eigLogger::log->trace("Starting eig_zgeev (lapacke)");
    using Scalar = std::complex<double>;
    // int lwork   =  2*2*L;
    // For some reason the recommended lwork from netlib doesn't work. It's better to ask lapack with a query.
    int lrwork  =  2*L;
    int info   = 0;
    Scalar lwork_query;
    std::vector<double> rwork  ( (unsigned long) lrwork);

    info = LAPACKE_zgeev_work(LAPACK_COL_MAJOR,'V','V',L,
                              reinterpret_cast< ComplexDouble_ *>(const_cast<Scalar *>(matrix)),
                              L,
                              reinterpret_cast< ComplexDouble_ *>(eigvals),
                              reinterpret_cast< ComplexDouble_ *>(eigvecsL),
                              L,
                              reinterpret_cast< ComplexDouble_ *>(eigvecsR),
                              L,
                              reinterpret_cast< ComplexDouble_ *>(&lwork_query),
                              -1,
                              rwork.data());
    int lwork = (int) std::real(2.0*lwork_query); //Make it twice as big for performance.
    std::vector<Scalar> work  ( (unsigned long)lwork );
    info = LAPACKE_zgeev_work(LAPACK_COL_MAJOR,'V','V',L,
                              reinterpret_cast< ComplexDouble_ *>(const_cast<Scalar *>(matrix)),
                              L,
                              reinterpret_cast< ComplexDouble_ *>(eigvals),
                              reinterpret_cast< ComplexDouble_ *>(eigvecsL),
                              L,
                              reinterpret_cast< ComplexDouble_ *>(eigvecsR),
                              L,
                              reinterpret_cast< ComplexDouble_ *>(work.data()),
                              lwork,
                              rwork.data());

    return info;
}








