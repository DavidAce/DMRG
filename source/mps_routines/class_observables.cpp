//
// Created by david on 2017-11-12.
//

#include <iomanip>
#include "class_observables.h"
#include "general/class_eig_spectra_wrapper.h"
#include "general/class_eig_arpack_wrapper.h"

class_observables::class_observables(class_superblock &superblockRef, SimulationType sim_)
                    :superblock(superblockRef),
                    sim(sim_)
{
    Sim2String = {
            {SimulationType::iDMRG,     "iDMRG: "},
            {SimulationType::fDMRG,     "fDMRG: "},
            {SimulationType::FES_iDMRG, "FES_iDMRG: "},
            {SimulationType::iTEBD,     "iTEBD: "},
            {SimulationType::FES_iTEBD, "FES_iTEBD: "}
    };
}


double class_observables::get_expectationvalue(const Tensor<double,4> &MPO){
    Tensor<double,4> theta  = superblock.MPS.get_theta();
    Tensor<double,0> result =
            superblock.Lblock.block
                    .contract(theta.conjugate(),idx<1>({0},{2}))
                    .contract(MPO,              idx<2>({1,2},{0,2}))
                    .contract(MPO,              idx<2>({3,1},{0,2}))
                    .contract(theta,            idx<3>({0,2,4},{2,0,1}))
                    .contract(superblock.Rblock.block,     idx<3>({0,1,2},{0,2,1}));
    return result(0)/ superblock.chain_length;
}

double class_observables::get_expectationvalue(const Tensor<std::complex<double>,4> &MPO){
    using T = std::complex<double>;
    Tensor<T,4> theta  = superblock.MPS.get_theta().cast<T>();
    Tensor<T,0> result =
            superblock.Lblock.block.cast<T>()
                    .contract(theta.conjugate(),idx<1>({0},{2}))
                    .contract(MPO,              idx<2>({1,2},{0,2}))
                    .contract(MPO,              idx<2>({3,1},{0,2}))
                    .contract(theta,            idx<3>({0,2,4},{2,0,1}))
                    .contract(superblock.Rblock.block.cast<T>(),     idx<3>({0,1,2},{0,2,1}));
    return result(0).real()/ superblock.chain_length;
}


double class_observables::get_second_cumulant(){

    array4 shape4 = superblock.shape4;
    long sizeL_MPO = shape4[1] * shape4[1] * superblock.H.F.dimension(0);
    long sizeR_MPO = shape4[3] * shape4[3] * superblock.H.F.dimension(1);
    long sizeL     = shape4[1] * shape4[1];
    long sizeR     = shape4[3] * shape4[3];
//
//    Tensor<std::complex<double>,2> transf_mat_R_MPO = superblock.MPS.get_transfer_2_site_matrix_R(superblock.H.GG).reshape(array2{sizeL_MPO,sizeR_MPO});
//    Tensor<std::complex<double>,2> transf_mat_R     = superblock.MPS.get_transfer_2_site_matrix_R(superblock.H.IGG).reshape(array2{sizeL_MPO    ,sizeL_MPO});

//    Tensor<double,2> transf_mat_R_MPO = superblock.MPS.get_transfer_matrix_R(superblock.H.F) .reshape(array2{sizeL_MPO,sizeR_MPO});
//    Tensor<double,2> transf_mat_R     = superblock.MPS.get_transfer_matrix_R(superblock.H.IF).reshape(array2{sizeL_MPO,sizeL_MPO});
//


//    class_eig_arpack eig;
//    auto largest_eigval_MPO = eig.solve_dominant<arpack::Form::SYMMETRIC, arpack::Ritz::LM, arpack::Side::R, false>(transf_mat_R_MPO, 1);
//    auto largest_eigval     = eig.solve_dominant<arpack::Form::SYMMETRIC, arpack::Ritz::LM, arpack::Side::R, false>(transf_mat_R    , 1);
////    cout << "\n hello " << endl;
//    cout << "\n lm = " << largest_eigval_MPO << endl;
//    cout << " l  = "   << largest_eigval     << endl;
//    double G_inf = std::pow(std::fabs(std::norm(largest_eigval_MPO(0))/std::norm(largest_eigval(0))), 0.5);
//    double G_inf = std::pow(std::fabs(std::norm(largest_eigval_MPO(0))/std::norm(largest_eigval(0))), 1.0);
    return 1.0;

}


//double class_observables::get_expectationvalue(const Tensor<Scalar,4> &MPO){
//    Tensor<Scalar,4> theta  = superblock.MPS.get_theta();
//    Tensor<Scalar,0> result =
//            superblock.Lblock.block
//            .contract(theta.conjugate(),idx<1>({0},{2}))
//            .contract(MPO,              idx<2>({1,2},{0,2}))
//            .contract(MPO,              idx<2>({3,1},{0,2}))
//            .contract(theta,            idx<3>({0,2,4},{2,0,1}))
//            .contract(superblock.Rblock.block,     idx<3>({0,1,2},{0,2,1}));
//    return static_cast<double>(result(0))/ superblock.chain_length;
//}


double class_observables::get_energy(){
    if (Sim2String[sim].find("TEBD") != std::string::npos ){
        auto theta  = superblock.MPS.get_theta();
        Tensor<double,0> E1     = superblock.H.H_asTensor
                                           .contract(theta.conjugate(), idx<2>({0,1},{0,1}))
                                           .contract(theta,             idx<4>({0,1,2,3},{0,1,2,3}));
        return E1(0);
    }else {
        return get_expectationvalue(superblock.H.M);
    }
}

double class_observables::get_entropy(){
    Tensor<Scalar,0> result1  = -superblock.MPS.LA.square()
            .contract(superblock.MPS.LA.square().log().eval(), idx<1>({0},{0}));
    return static_cast<double>(result1(0)) ;
}


//Tensor<std::complex<double>,4> class_observables::get_Lblock_sq(){
//    class_eig<Scalar> eig;
//    using T = std::complex<double>;
//    if (superblock.MPS.LA.size() == superblock.MPS.LB.size()) {
//        array4 shape = {
//                superblock.MPS.GA.dimension(2),
//                superblock.MPS.GA.dimension(2),
//                superblock.H.M.dimension(1),
//                superblock.H.M.dimension(1)
//        };
//        long sizeL = superblock.MPS.GA.dimension(1)
//                     * superblock.H.M.dimension(0)
//                     * superblock.H.M.dimension(0)
//                     * superblock.MPS.GA.dimension(1);
//        long sizeR = superblock.MPS.GA.dimension(2)
//                     * superblock.H.M.dimension(1)
//                     * superblock.H.M.dimension(1)
//                     * superblock.MPS.GA.dimension(2);
//        Tensor<Scalar, 2> LGMPOGL = asDiagonal(superblock.MPS.L_tail).cast<T>()
//                .contract(superblock.MPS.GA.cast<T>(), idx<1>({1}, {1}))
//                .contract(superblock.H.M, idx<1>({1}, {2}))
//                .contract(superblock.H.M, idx<1>({4}, {2}))
//                .contract(superblock.MPS.GA.cast<T>(), idx<1>({6}, {0}))
//                .contract(asDiagonal(superblock.MPS.L_tail).cast<T>(), idx<1>({6}, {1}))
//                .shuffle(array8{0, 7, 2, 4, 1, 6, 3, 5})
//                .reshape(array2{sizeL, sizeR}).real();
//        auto[H2, H2val] = eig.solve_dominant<eig_Form::GENERAL, eig_Mode::LARGEST_MAGN, eig_Side::R>(LGMPOGL, 1,{sizeR, 1});
//        return H2.real().reshape(shape);
//    }
//};
//
//Tensor<std::complex<double>,4> class_observables::get_Rblock_sq(){
//    class_eig<Scalar> eig;
//    using T = std::complex<double>;
//    if (superblock.MPS.LA.size() == superblock.MPS.LB.size()) {
//        array4 shape = {
//                superblock.MPS.GA.dimension(1),
//                superblock.MPS.GA.dimension(1),
//                superblock.H.M.dimension(0),
//                superblock.H.M.dimension(0)
//        };
//        long sizeL = superblock.MPS.GB.dimension(1)
//                     * superblock.H.M.dimension(0)
//                     * superblock.H.M.dimension(0)
//                     * superblock.MPS.GB.dimension(1);
//        long sizeR = superblock.MPS.GB.dimension(2)
//                     * superblock.H.M.dimension(1)
//                     * superblock.H.M.dimension(1)
//                     * superblock.MPS.GB.dimension(2);
//        Tensor<Scalar, 2> LGMPOGL = asDiagonal(superblock.MPS.LB).cast<T>()
//                .contract(superblock.MPS.GB.cast<T>(), idx<1>({0}, {2}))
//                .contract(superblock.H.M, idx<1>({1}, {2}))
//                .contract(superblock.H.M, idx<1>({4}, {2}))
//                .contract(superblock.MPS.GA.cast<T>(), idx<1>({6}, {0}))
//                .contract(asDiagonal(superblock.MPS.LB).cast<T>(), idx<1>({7}, {0}))
//                .shuffle(array8{0, 7, 3, 5, 1, 6, 3, 5})
//                .reshape(array2{sizeL, sizeR}).real();
//
//        auto[H2, H2val] = eig.solve_dominant<eig_Form::GENERAL, eig_Mode::LARGEST_MAGN, eig_Side::L>(LGMPOGL, 1,{sizeL, 1});
//        return H2.reshape(shape);
//    }
//};

double class_observables::get_variance(){


    if (Sim2String[sim].find("TEBD") != std::string::npos ){
        Tensor<double,4> theta = superblock.MPS.get_theta();
        Tensor<double,4> H_sq  = superblock.H.H_asTensor.contract(superblock.H.H_asTensor, idx<2>({0,1},{2,3}));
        Tensor<double,0> Var   =  H_sq.contract(theta.conjugate(), idx<2>({0,1},{0,1}))
                .contract(theta,             idx<4>({0,1,2,3},{0,1,2,3}));
        return ( Var(0)- std::pow(get_energy(), 2));
    }else{
        if (superblock.MPS.LA.size() == superblock.MPS.LB.size()) {
            return get_second_cumulant();
        }else{
            return 1;
        }
    }


}

double class_observables::get_truncation_error(){
    return superblock.truncation_error;
}

long class_observables::get_chi(){
    return superblock.chi;
}

long class_observables::get_chi_max(){
    return superblock.chi_max;
}

long class_observables::get_chain_length(){
    return superblock.chain_length;
}

void class_observables::print_status_full(int verbosity) {
    if (verbosity >= 1) {
        cout << setprecision(16) << fixed << left;
        cout << setw(20) << "E_inf_" + Sim2String.at(sim) << " = "   << get_energy() << '\n';
        if(verbosity >= 2) {
            cout  << setw(20) << "Entanglement Entropy = " << setprecision(16) << fixed      << get_entropy() << '\n';
            cout  << setw(20) << "Truncation error     = " << setprecision(4)  << scientific << get_truncation_error() << '\n';
            cout  << setw(20) << "Variance             = " << setprecision(16) << scientific << get_variance() << '\n';
            cout  << setw(20) << "chi_max              = " << setprecision(4)  << fixed      << get_chi_max() << '\n';
            cout  << setw(20) << "chi                  = " << setprecision(4)  << fixed      << get_chi()     << '\n';

        }
    }
}


void class_observables::print_status_update(int step) {
    cout << left  << Sim2String.at(sim)
         << left  << setw(10) << step
         << left  << "E: "                      << setw(25) << setprecision(16)    << fixed   << get_energy()
         << left  << "Var(E): "                 << setw(25) << setprecision(16)    << fixed   << get_variance()
         << left  << "Entanglement Entropy: "   << setw(25) << setprecision(16)    << fixed   << get_entropy()
         << left  << "SVD truncation: "         << setw(25) << setprecision(16)    << fixed   << get_truncation_error()
         << left  << "\u03C7_max: "             << setw(3)  << setprecision(4)     << fixed   << get_chi_max()
         << left  << "\u03C7: "                 << setw(20) << setprecision(4)     << fixed   << get_chi()
         << '\n';
}
