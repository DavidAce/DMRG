
#include <complex>
//
#ifdef MKL_AVAILABLE
#include <mkl_service.h>
//#define MKL_Complex8 std::complex<float>
//#define MKL_Complex16 std::complex<double>
#endif
//
////#include <mkl.h>
////#define LAPACKF_H
////#define LAPACKC_H
#include <arpack++/arcomp.h>
#include <arpack++/ardscomp.h>
#include <arpack++/ardnsmat.h>




int main(int argc, char* argv[])
{

#ifdef MKL_VERBOSE
    MKL_Verbose(1);
#endif

    int                n = 4;
    int                nev = 2;
    int                ncv = 4;
    char ritz[10];
    std::string("LR").copy(ritz,2);
    using Scalar = std::complex<double>;
//    using Scalar = arcomplex<double>;

    std::vector<Scalar> data = {Scalar(1.0,+0.0),Scalar(0.0,+1.0),Scalar(0.0,+1.0),Scalar(0.0,+1.0),
                                Scalar(0.0,-1.0),Scalar(2.0,+0.0),Scalar(0.0,+2.0),Scalar(0.0,+2.0),
                                Scalar(0.0,-1.0),Scalar(0.0,+2.0),Scalar(3.0,+3.0),Scalar(0.0,+3.0),
                                Scalar(0.0,-1.0),Scalar(0.0,+2.0),Scalar(0.0,-3.0),Scalar(4.0,+0.0)};

    // Creating a complex matrix.
    ARdsNonSymMatrix<Scalar, double> A(n, data.data());

    ARluCompStdEig<double> solver(nev, A,ritz,ncv);
    std::cout << "HERE " << std::endl;

    solver.FindEigenvectors();
    std::cout << "Eigenvalues: " << solver.Eigenvalue(0)  << " " << solver.Eigenvalue(1) << std::endl;
    std::cout << "Iterations : " << solver.GetIter()   << std::endl;
    if (solver.EigenvaluesFound()){return 0;}
    else {return 1;}
}