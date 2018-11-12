//
// Created by david on 2018-10-29.
//

#include "class_eigsolver.h"


//using namespace eigSetting;




void class_eigsolver::conf_init  (const int L,
                                  const int nev,
                                  const int ncv,
                                  const std::complex<double> sigma          ,
                                  const eigutils::eigSetting::Type type      ,
                                  const eigutils::eigSetting::Form form      ,
                                  const eigutils::eigSetting::Ritz ritz      ,
                                  const eigutils::eigSetting::Side side      ,
                                  const eigutils::eigSetting::Storage storage,
                                  const bool compute_eigvecs_               ,
                                  const bool remove_phase_
)
{
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


















