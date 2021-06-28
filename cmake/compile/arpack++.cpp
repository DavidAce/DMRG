#include <complex>
#include <vector>
#if __has_include(<arpackpp/arcomp.h>)
    #include <arpackpp/arcomp.h>
    #include <arpackpp/ardscomp.h>
    #include <arpackpp/ardnsmat.h>
#elif __has_include(<arpack++/arcomp.h>)
    #include <arpack++/arcomp.h>
    #include <arpack++/ardscomp.h>
    #include <arpack++/ardnsmat.h>
#elif __has_include(<arcomp.h>)
    #include <arcomp.h>
    #include <ardscomp.h>
    #include <ardnsmat.h>
#else
#error Could not include arpack headers correctly
#endif
int main(){

    int                n = 4;
    int                nev = 2;
    int                ncv = 4;
    char ritz[10];
    std::string("LR").copy(ritz,2);
    using Scalar = std::complex<double>;
    std::vector<Scalar> data = {Scalar(1.0,+0.0),Scalar(0.0,+1.0),Scalar(0.0,+1.0),Scalar(0.0,+1.0),
                                Scalar(0.0,-1.0),Scalar(2.0,+0.0),Scalar(0.0,+2.0),Scalar(0.0,+2.0),
                                Scalar(0.0,-1.0),Scalar(0.0,+2.0),Scalar(3.0,+3.0),Scalar(0.0,+3.0),
                                Scalar(0.0,-1.0),Scalar(0.0,+2.0),Scalar(0.0,-3.0),Scalar(4.0,+0.0)};
    ARdsNonSymMatrix<Scalar, double> A(n, data.data());
    ARluCompStdEig<double> solver(nev, A,ritz,ncv);
    solver.FindEigenvectors();
    if (solver.EigenvaluesFound()){return 0;}
    else {return 1;}
}