//
// Created by david on 2018-10-17.
//

#ifndef DMRG_CLASS_OPTIMIZATION_CPPOPTLIB_H
#define DMRG_CLASS_OPTIMIZATION_CPPOPTLIB_H

#include <Eigen/Core>
#include <cppoptlib/meta.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/solver/bfgssolver.h>
#include <cppoptlib/solver/lbfgsbsolver.h>
#include <cppoptlib/solver/conjugatedgradientdescentsolver.h>
#include <cppoptlib/solver/newtondescentsolver.h>
#include <cppoptlib/solver/gradientdescentsolver.h>
#include <cppoptlib/solver/neldermeadsolver.h>
#include <cppoptlib/solver/cmaessolver.h>


enum class opt_algorithm {BFGS,LBFGS, GD,CGD, ND, NM, CM};



class class_optimization_cppoptlib {
public:
    cppoptlib::Criteria<double> criteria_obj = cppoptlib::Criteria<double>::defaults();
    cppoptlib::Status status_obj = cppoptlib::Status::NotStarted;
    class_optimization_cppoptlib(){
        criteria_obj.iterations = 1000;                                                     // Increase the number of allowed iterations
        criteria_obj.fDelta     = 0;
        criteria_obj.xDelta     = 0;
        criteria_obj.gradNorm   = 0;
        criteria_obj.condition  = 0;
    };



    template <opt_algorithm Algo = opt_algorithm::BFGS,typename FunctorType, typename FVectorType>
    void minimize(FunctorType &objFun,FVectorType &x0){
        if constexpr(Algo == opt_algorithm::BFGS){
            cppoptlib::BfgsSolver<FunctorType> solver;
            solver.setStopCriteria(criteria_obj);
            solver.minimize(objFun, x0);
            status_obj   = solver.status();
            criteria_obj = solver.criteria();
        }
        else
        if constexpr(Algo == opt_algorithm::LBFGS){
            cppoptlib::LbfgsbSolver<FunctorType> solver;
            solver.setStopCriteria(criteria_obj);
            solver.minimize(objFun, x0);
            status_obj   = solver.status();
            criteria_obj = solver.criteria();
        }
        else
        if constexpr(Algo == opt_algorithm::GD){
            cppoptlib::GradientDescentSolver<FunctorType> solver;
            solver.setStopCriteria(criteria_obj);
            solver.minimize(objFun, x0);
            status_obj   = solver.status();
            criteria_obj = solver.criteria();
        }
        else
        if constexpr(Algo == opt_algorithm::CGD){
            cppoptlib::ConjugatedGradientDescentSolver<FunctorType> solver;
            solver.setStopCriteria(criteria_obj);
            solver.minimize(objFun, x0);
            status_obj   = solver.status();
            criteria_obj = solver.criteria();
        }
        else
        if constexpr(Algo == opt_algorithm::ND){
            cppoptlib::NewtonDescentSolver<FunctorType> solver;
            solver.setStopCriteria(criteria_obj);
            solver.minimize(objFun, x0);
            status_obj   = solver.status();
            criteria_obj = solver.criteria();
        }
        else
        if constexpr(Algo == opt_algorithm::NM){
            cppoptlib::NelderMeadSolver<FunctorType> solver;
            solver.setStopCriteria(criteria_obj);
            solver.minimize(objFun, x0);
            status_obj   = solver.status();
            criteria_obj = solver.criteria();
        }
        else
        if constexpr(Algo == opt_algorithm::CM){
            cppoptlib::CMAesSolver<FunctorType> solver;
            solver.setStopCriteria(criteria_obj);
            solver.minimize(objFun, x0);
            status_obj   = solver.status();
            criteria_obj = solver.criteria();
        }
    }

    const auto &status(){return status_obj;}
    const auto &criteria(){return criteria_obj;}
};


#endif //DMRG_CLASS_OPTIMIZATION_CPPOPTLIB_H
