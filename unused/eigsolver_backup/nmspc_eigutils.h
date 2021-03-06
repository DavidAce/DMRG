//
// Created by david on 2018-06-07.
//

#pragma once


#include <iostream>
#include <array>
#include <map>
#include <complex>
#include <vector>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace eigutils{

    namespace eigSetting{
        enum class Form{SYMMETRIC, NONSYMMETRIC};       // Real Symmetric, Real General or Complex General
        enum class Storage {DENSE,SPARSE,STL};          // Eigen Dense or sparse, or std::vector for container
        enum class Shift {ON,OFF};                      // Enable or disable shift invert
        enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE};   // Choice of eigenvalue. LA is largest algebraic, and so on.
        enum class Side {L,R};                          // Left or right eigenvectors
        enum class Type {REAL,CPLX};                    // Real or complex, i.e. double or std::complex<double> matrix

        inline Ritz stringToRitz(std::string_view ritzstring){
            if (ritzstring == "LA") return Ritz::LA;
            if (ritzstring == "SA") return Ritz::SA;
            if (ritzstring == "LM") return Ritz::LM;
            if (ritzstring == "SM") return Ritz::SM;
            if (ritzstring == "LR") return Ritz::LR;
            if (ritzstring == "SR") return Ritz::SR;
            if (ritzstring == "LI") return Ritz::LI;
            if (ritzstring == "SI") return Ritz::SI;
            if (ritzstring == "BE") return Ritz::BE;
            throw std::runtime_error("Wrong ritz string: " + std::string(ritzstring));
        }

        inline std::string_view RitzToString(Ritz ritz){
            if (ritz == Ritz::LA) return  "LA";
            if (ritz == Ritz::SA) return  "SA";
            if (ritz == Ritz::LM) return  "LM";
            if (ritz == Ritz::SM) return  "SM";
            if (ritz == Ritz::LR) return  "LR";
            if (ritz == Ritz::SR) return  "SR";
            if (ritz == Ritz::LI) return  "LI";
            if (ritz == Ritz::SI) return  "SI";
            if (ritz == Ritz::BE) return  "BE";
            throw std::runtime_error("Wrong ritz enum");
        }


    }

    class eigConfig{
    public:
        bool confOK = false;
        eigSetting::Form           form        = eigSetting::Form::NONSYMMETRIC;
        eigSetting::Storage        storage     = eigSetting::Storage::DENSE;
        eigSetting::Shift          shift       = eigSetting::Shift::OFF;
        eigSetting::Side           side        = eigSetting::Side::R;
        eigSetting::Ritz           ritz        = eigSetting::Ritz::LM;
        eigSetting::Type           type        = eigSetting::Type::REAL;
        bool    compute_eigvecs                = true;
        bool    remove_phase                   = false;
        double  eigThreshold                   = 1e-12;
        int     eigMaxIter                     = 2000;
        int     eigMaxNev                      = 1;
        int     eigMaxNcv                      = 16;
        std::complex<double>  sigma            = std::numeric_limits<std::complex<double>>::quiet_NaN();     // Sigma value for shift-invert mode.
        std::string get_ritz(){return std::string(eigutils::eigSetting::RitzToString(ritz));}
        void checkRitz()
        // Checks that ritz is valid for the given problem.
        // The valid ritzes are stated in the arpack++ manual page 78.
        {
            using namespace eigSetting;
            std::string suggestion;
            if (type==Type::CPLX or form==Form::NONSYMMETRIC){
                if (ritz==Ritz::LA or
                    ritz==Ritz::SA or
                    ritz==Ritz::BE)
                {
                    std::cerr << "WARNING: Invalid ritz for cplx or nonsym problem: " << RitzToString(ritz) << std::endl;
                    if (ritz==Ritz::LA) suggestion = "LR";
                    if (ritz==Ritz::SA) suggestion = "SR";
                    if (ritz==Ritz::BE) suggestion = "LM";
                    std::cerr << "         Suggested ritz: " << suggestion << std::endl;
                }
            }else if (type==Type::REAL and form==Form::SYMMETRIC) {
                if (ritz==Ritz::LR or
                    ritz==Ritz::SR or
                    ritz==Ritz::LI or
                    ritz==Ritz::SI)
                {
                    std::cerr << "WARNING: Invalid ritz for real and sym problem: " << RitzToString(ritz)<< std::endl;
                    if (ritz==Ritz::LR) suggestion = "LA";
                    if (ritz==Ritz::SR) suggestion = "SA";
                    if (ritz==Ritz::LI) suggestion = "LM";
                    if (ritz==Ritz::SI) suggestion = "SM";
                    std::cerr << "         Suggested ritz: " << suggestion << std::endl;
                }
            }
            if(not suggestion.empty())
                throw std::runtime_error("Invalid ritz");
            confOK = true;
        }

    };


    class eigSolution{
    public:
        using Scalar = std::complex<double>;
        // For symmetric problems
        std::vector<double> eigvecs_real;
        std::vector<Scalar> eigvecs_cplx;
        std::vector<double> eigvals_real;

        // For nonsymmetric problems
        std::vector<double> eigvecsR_real;
        std::vector<double> eigvecsL_real;
        std::vector<Scalar> eigvecsR_cplx;
        std::vector<Scalar> eigvecsL_cplx;
        std::vector<Scalar> eigvals_cplx;


        template<eigutils::eigSetting::Type   type,
                eigutils::eigSetting::Form   form,
                eigutils::eigSetting::Side   side = eigutils::eigSetting::Side::R >
        auto & get_eigvecs(){
            using namespace eigutils::eigSetting;
            if constexpr(form == Form::SYMMETRIC){
                if constexpr(type == Type::REAL)                    {return eigvecs_real;}
                if constexpr(type == Type::CPLX)                    {return eigvecs_cplx;}
            }else if constexpr(form == Form::NONSYMMETRIC and side == Side::R){
                if constexpr(type == Type::REAL and side == Side::R){return eigvecsR_cplx;}
                if constexpr(type == Type::CPLX and side == Side::R){return eigvecsR_cplx;}
            }else if constexpr(form == Form::NONSYMMETRIC and side == Side::L){
                if constexpr(type == Type::REAL and side == Side::L){return eigvecsL_cplx;}
                if constexpr(type == Type::CPLX and side == Side::L){return eigvecsL_cplx;}
            }
        }

        template<typename Scalar,
                eigutils::eigSetting::Form   form,
                eigutils::eigSetting::Side   side = eigutils::eigSetting::Side::R >
        auto & get_eigvecs(){
            using namespace eigutils::eigSetting;
            if constexpr (std::is_same<double, Scalar>::value){
                return get_eigvecs<Type::REAL,form,side>();
            }else if constexpr (std::is_same<std::complex<double>, Scalar>::value){
                return get_eigvecs<Type::CPLX,form,side>();
            }
        }


        template<eigutils::eigSetting::Form form>
        auto & get_eigvals(){
            using namespace eigutils::eigSetting;
            if constexpr(form == Form::SYMMETRIC)   {return eigvals_real;}
            if constexpr(form == Form::NONSYMMETRIC){return eigvals_cplx;}
        }


        struct Meta{
            int     rows            = 0;
            int     cols            = 0;
            int     iter            = 0;
            int     nev             = 0; // Found eigenvectors. aka cols.
            int     nev_converged   = 0;
            int     n               = 0; // Linear dimension of the input matrix to diagonalize, aka rows.
            int     counter         = 0;
            int     ncv_used        = 0;
            bool    eigvals_found   = false; // For all problems
            bool    eigvecs_found   = false; // For symmetric problems
            bool    eigvecsR_found  = false; // For nonsymmetric problems
            bool    eigvecsL_found  = false; // For nonsymmetric problems
            eigSetting::Form  form;
            eigSetting::Ritz  ritz;
            eigSetting::Type  type;
            eigSetting::Side  side;
        } meta;
        void reset(){
            eigvals_real.clear();
            eigvals_cplx.clear();
            eigvecs_real.clear();
            eigvecs_cplx.clear();
            eigvecsR_real.clear();
            eigvecsL_real.clear();
            eigvecsR_cplx.clear();
            eigvecsL_cplx.clear();
            meta = Meta();
        }
    };



    namespace eigLogger{

        inline std::shared_ptr<spdlog::logger> log;
        inline void enableTimeStamp(std::shared_ptr<spdlog::logger> &log){
            if(log != nullptr) {
                log->set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
            }
        }
        inline void disableTimeStamp(std::shared_ptr<spdlog::logger> &log){
            if(log != nullptr){
                log->set_pattern("[%n]%^[%=8l]%$ %v");
            }
        }

        inline void setLogLevel(size_t levelZeroToSix){
            if (levelZeroToSix > 6) {
                throw std::runtime_error("Expected verbosity level integer in [0-6]. Got: " + std::to_string(levelZeroToSix));
            }
            auto lvlEnum = static_cast<spdlog::level::level_enum>(levelZeroToSix);

            // Set console settings
            log->set_level(lvlEnum);
        }

        inline void setLogger(std::string name, size_t levelZeroToSix = 2, bool timestamp = true){
            if(spdlog::get(name) == nullptr){
                log = spdlog::stdout_color_mt(name);
                if (timestamp){enableTimeStamp(log);}
                else{disableTimeStamp(log); }
                setLogLevel(levelZeroToSix);
            }else{
                log = spdlog::get(name);
            }
        }

    }

}
