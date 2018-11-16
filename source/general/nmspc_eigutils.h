//
// Created by david on 2018-06-07.
//

#ifndef ARPACK_EXTRA_H
#define ARPACK_EXTRA_H

#include <iostream>
#include <array>
#include <map>
#include <complex>
#include <vector>


namespace eigutils{

    namespace eigSetting{
        enum class Form{SYMMETRIC, NONSYMMETRIC};       // Real Symmetric, Real General or Complex General
        enum class Storage {DENSE,SPARSE,STL};          // Eigen Dense or sparse, or std::vector for container
        enum class Shift {ON,OFF};                      // Enable or disable shift invert
        enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE};   // Choice of eigenvalue. LA is largest algebraic, and so on.
        enum class Side {L,R};                          // Left or right eigenvectors
        enum class Type {REAL,CPLX};                    // Real or complex, i.e. double or std::complex<double> matrix
    }

    class eigConfig{
    public:
        bool confOK = false;
        using  MapType = std::map<eigSetting::Ritz, std::string>;
        MapType RitzToString;
        char ritz_char[3];


        eigSetting::Form           form        = eigSetting::Form::NONSYMMETRIC;
        eigSetting::Storage        storage     = eigSetting::Storage::DENSE;
        eigSetting::Shift          shift       = eigSetting::Shift::OFF;
        eigSetting::Side           side        = eigSetting::Side::R;
        eigSetting::Ritz           ritz        = eigSetting::Ritz::LM;
        eigSetting::Type           type        = eigSetting::Type::REAL;
        bool    compute_eigvecs                = false;
        bool    remove_phase                   = false;
        double  eigThreshold                   = 1e-12;
        int     eigMaxIter                     = 2000;
        int     eigMaxNev                      = 1;
        int     eigMaxNcv                      = 16;
        std::complex<double>  sigma            = std::numeric_limits<std::complex<double>>::quiet_NaN();     // Sigma value for shift-invert mode.

        eigConfig() {
            RitzToString = {
                    {eigSetting::Ritz::LA, "LA"},
                    {eigSetting::Ritz::SA, "SA"},
                    {eigSetting::Ritz::LM, "LM"},
                    {eigSetting::Ritz::SM, "SM"},
                    {eigSetting::Ritz::LR, "LR"},
                    {eigSetting::Ritz::SR, "SR"},
                    {eigSetting::Ritz::LI, "LI"},
                    {eigSetting::Ritz::SI, "SI"},
                    {eigSetting::Ritz::BE, "BE"}
            };
        }

        void writeRitzChar()
        // Writes ritz to string and checks that it is valid for the given problem.
        // The valid ritzes are stated in the arpack++ manual page 78.
        {
            using namespace eigSetting;
            if (type==Type::CPLX or form==Form::NONSYMMETRIC){
                if (ritz==Ritz::LA or
                    ritz==Ritz::SA or
                    ritz==Ritz::BE
                    )
                {
                    std::cerr << "WARNING: Invalid ritz for nonsym problem: " << RitzToString.at(ritz) << std::endl;
                    if (ritz==Ritz::LA){ritz = Ritz::LR;}
                    if (ritz==Ritz::SA){ritz = Ritz::SR;}
                    if (ritz==Ritz::BE){ritz = Ritz::LM;}
                    std::cerr << "         Changed ritz to : " << RitzToString.at(ritz)<< std::endl;
                }
            }else if (type==Type::REAL and form==Form::SYMMETRIC) {
                if (ritz==Ritz::LR or
                    ritz==Ritz::SR or
                    ritz==Ritz::LI or
                    ritz==Ritz::SI
                    )
                {
                    std::cerr << "WARNING: Invalid ritz for nonsym problem: " << RitzToString.at(ritz)<< std::endl;
                    if (ritz==Ritz::LR){ritz = Ritz::LA;}
                    if (ritz==Ritz::SR){ritz = Ritz::SA;}
                    if (ritz==Ritz::LI){ritz = Ritz::LM;}
                    if (ritz==Ritz::SI){ritz = Ritz::SM;}
                    std::cerr << "         Changed ritz to : " << RitzToString.at(ritz)<< std::endl;
                }
            }

            RitzToString.at(ritz).copy(ritz_char, 3);
            confOK = true;
        }

    };


    class eigSolution{
    public:
        using Scalar = std::complex<double>;
        std::vector<Scalar> eigvecs;
        std::vector<Scalar> eigvals;

        struct Meta{
            int     rows            = 0;
            int     cols            = 0;
            int     iter            = 0;
            int     nev             = 0; // Found eigenvectors. aka cols.
            int     nev_converged   = 0;
            int     n               = 0; // Linear dimension of the input matrix to diagonalize, aka rows.
            int     counter         = 0;
            int     ncv_used        = 0;
            bool    eigvecs_found   = false;
            bool    eigvals_found   = false;
        } meta;
        void reset(){
            eigvecs.clear();
            eigvals.clear();
            meta = Meta();
        }
    };

}



#endif //ARPACK_EXTRA_H
