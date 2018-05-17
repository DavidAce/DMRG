#include <arpack++/arcomp.h>
#include <arpack++/ardscomp.h>
#include <arpack++/ardnsmat.h>


int main()
{

    int                n = 4;
    int                nev = 2;
    int                ncv = 4;
    std::vector<std::complex<double>> data = {1.0+0.0i,0.0+1.0i,0.0+1.0i,0.0+1.0i,
                                              0.0-1.0i,2.0+0.0i,0.0+2.0i,0.0+2.0i,
                                              0.0-1.0i,0.0+2.0i,3.0+3.0i,0.0+3.0i,
                                              0.0-1.0i,0.0+2.0i,0.0-3.0i,4.0+0.0i};

    // Creating a complex matrix.
    ARdsNonSymMatrix<std::complex<double>, double> A(n, data.data());

    ARluCompStdEig<double> solver(nev, A,"LM",ncv);
    solver.FindEigenvectors();
    std::cout << "Eigenvalues: " << solver.Eigenvalue(0)  << " " << solver.Eigenvalue(1) << std::endl;
    std::cout << "Iterations : " << solver.GetIter()   << std::endl;
    if (solver.EigenvaluesFound()){return 0;}
    else {return 1;}
}