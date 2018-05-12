//
// Created by david on 2017-11-12.
//

#ifndef DMRG_CLASS_OBSERVABLES_H
#define DMRG_CLASS_OBSERVABLES_H
#include <map>
#include <memory>
#include <complex>
#include <IO/class_custom_cout.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/class_tic_toc.h>

using namespace std::complex_literals;

class class_superblock;
class class_mps;
class class_finite_chain_storage;
/*!
 * \class class_measurement
 * \brief A class for measuring observables
 * This class extracts observables, expectation values, like energy and entropy from the MPS, as well as other useful numbers like \f$\chi\f$ and truncation errors.
*/

class class_measurement {
public:
    using Scalar = std::complex<double>;
private:
    std::shared_ptr<class_superblock> superblock;
    std::shared_ptr<class_finite_chain_storage> env_storage;
    class_custom_cout ccout;
    Scalar moment_generating_function(std::shared_ptr<class_mps> MPS_original,
                                                       std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec);

//    Scalar moment_generating_function_2(std::shared_ptr<class_mps> MPS_original,
//                                      std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec);
//

//    class_parity_mpo parity_mpo;
    Scalar a0 = (0.0 + 1.0i) *0.0 ;
    Scalar a  = (0.0 + 1.0i) *5e-2;
    Scalar b  = (0.0 + 1.0i) *1e-2;
    Scalar c  = (0.0 + 1.0i) *5e-3;
    Scalar d  = (0.0 + 1.0i) *1e-3;
    std::vector<Eigen::Tensor<Scalar,4>> mom_vec0;
    std::vector<Eigen::Tensor<Scalar,4>> mom_vecA;
    std::vector<Eigen::Tensor<Scalar,4>> mom_vecB;
    std::vector<Eigen::Tensor<Scalar,4>> mom_vecC;
    std::vector<Eigen::Tensor<Scalar,4>> mom_vecD;

private:
    double compute_energy_MPO();
    double compute_energy_H();
    double compute_entanglement_entropy();
    double compute_infinite_variance_MPO();
    double compute_infinite_variance_H();
    double compute_parity();
    std::pair<double,double> compute_infinite_moments_G(Scalar a, std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec);


    double energy1   = 0; double energy2   = 0; double energy3   = 0; double energy4   = 0; double energy5   = 0; double energy6   = 0;
    double variance1 = 0; double variance2 = 0; double variance3 = 0; double variance4 = 0; double variance5 = 0; double variance6 = 0;
    double entanglement_entropy = 0;
    double truncation_error = 0;
    double parity = 0;
    Eigen::Tensor<Scalar,0> E_evn, E_odd;

    public:
    SimulationType sim;
    bool is_measured = false;

    explicit class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_);
    void   compute_all_observables_from_superblock();
    void   compute_all_observables_from_finite_chain();

    double get_energy1();               /*! Computes the current energy.*/
    double get_energy2();               /*! Computes the current energy.*/
    double get_energy3();               /*! Computes the current energy.*/
    double get_energy4();               /*! Computes the current energy.*/
    double get_energy5();               /*! Computes the current energy.*/
    double get_energy6();               /*! Computes the current energy.*/


    double get_variance1();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance2();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance3();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance4();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance5();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance6();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_entanglement_entropy();
    double get_truncation_error();
    double get_parity();
    long   get_chi();
    long   get_chain_length();


    // Profiling
    class_tic_toc t_ene_mpo;
    class_tic_toc t_ene_ham;
    class_tic_toc t_ene_gen;
    class_tic_toc t_var_mpo;
    class_tic_toc t_var_ham;
    class_tic_toc t_var_gen;
    class_tic_toc t_entropy;
    class_tic_toc t_temp1;
    class_tic_toc t_temp2;
    class_tic_toc t_temp3;
    class_tic_toc t_temp4;

    void set_profiling_labels();
    void print_profiling(class_tic_toc &t_parent);
};


#endif //DMRG_CLASS_OBSERVABLES_H
