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
class class_mps_2site;
class class_mps_util;
class class_finite_chain_sweeper;
/*!
 * \class class_measurement
 * \brief A class for measuring observables
 * This class extracts observables, expectation values, like energy and entropy from the MPS, as well as other useful numbers like \f$\chi\f$ and truncation errors.
*/

class class_measurement {
public:
    using Scalar = std::complex<double>;
private:
    std::shared_ptr<const class_superblock>           superblock;
    SimulationType sim_type;
    std::shared_ptr<const class_finite_chain_sweeper> env_storage;
    std::shared_ptr<class_mps_util> mps_util;

    class_custom_cout ccout;
    Scalar moment_generating_function(const std::unique_ptr<class_mps_2site> &MPS_original,
                                                       std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec);

//    class_parity_mpo parity_mpo;
    Scalar a  = (0.0 + 1.0i) *5e-3;
    std::vector<Eigen::Tensor<Scalar,4>> mom_vecA;

    Eigen::MatrixXcd h_evn;
    Eigen::MatrixXcd h_odd;


    void compute_energy_mpo();
    void compute_energy_ham();
    void compute_entanglement_entropy();
    void compute_energy_variance_mpo();
    void compute_energy_variance_ham();
    void compute_energy_and_variance_mom(Scalar a, std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec);

    void compute_parity();


    double energy_mpo_all_sites;
    double variance_mpo_all_sites;

    double energy_mpo   = 0; double energy_ham   = 0; double energy_mom   = 0; //double energy4   = 0; double energy5   = 0; double energy6   = 0;
    double variance_mpo = 0; double variance_ham = 0; double variance_mom = 0; //double variance4 = 0; double variance5 = 0; double variance6 = 0;

    double norm_chain   = 0;
    double energy_chain = 0;
    double variance_chain = 0;

    double entanglement_entropy = 0;
    double parity = 0;
    Eigen::Tensor<Scalar,0> E_evn, E_odd;
    bool is_measured = false;

    public:
    explicit class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_);
    explicit class_measurement(std::shared_ptr<class_superblock> superblock_, std::shared_ptr<class_finite_chain_sweeper> env_storage_, SimulationType sim_);
    void   compute_all_observables_from_superblock();
    void   compute_all_observables_from_superblock(const Eigen::Tensor<Scalar,4> &theta);
    void   compute_all_observables_from_finite_chain();
    Eigen::Tensor<Scalar,1> mps_chain;

    double get_energy_mpo();               /*! Computes the current energy.*/
    double get_energy_ham();               /*! Computes the current energy.*/
    double get_energy_mom();               /*! Computes the current energy.*/
    double get_energy_mpo_all_sites();

    double get_variance_mpo();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance_ham();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance_mom();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
//  double get_variance_mpo_all_sites();
    double get_entanglement_entropy();
    double get_truncation_error();
    double get_parity();
    long   get_chi();
    long   get_chain_length();

    void compute_finite_chain_norm();
    void compute_finite_chain_norm2();
    void compute_finite_chain_energy();
    void compute_finite_chain_energy_variance();
    void compute_finite_chain_mps_state();


    void set_not_measured(){is_measured = false;}
    bool has_been_measured(){return is_measured;}
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
