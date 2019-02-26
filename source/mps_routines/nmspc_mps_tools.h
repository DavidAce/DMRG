//
// Created by david on 2019-01-29.
//

#ifndef NMSPC_FINITE_CHAIN_TOOLS_H
#define NMSPC_FINITE_CHAIN_TOOLS_H
#include <memory>
#include <string>
#include <general/nmspc_tensor_extra.h>
class class_superblock;
class class_mps_2site;
class class_finite_chain_state;
class class_hdf5_file;
class class_hamiltonian_base;
class class_simulation_state;

namespace MPS_Tools{

    namespace Finite
    /*!
     * Functions for finite MPS algorithms like fDMRG and xDMRG
     * These functions require the class_finite_chain_state containing the
     * MPS and environments for all sites of the lattice.
     */
    {
        namespace Chain {
            extern void initialize_state(class_finite_chain_state &state,std::string model_type, const size_t length, const size_t seed=0);
            extern void randomize_mpos  (class_finite_chain_state &state, const size_t seed);

            extern void copy_superblock_to_state(class_finite_chain_state &state, const class_superblock &superblock);    /*!< Update the MPS, MPO and ENV stored at current position.*/
            extern void copy_superblock_mps_to_state(class_finite_chain_state &state,const class_superblock &superblock);    /*!< Update the MPS stored at current position.*/
            extern void copy_superblock_mpo_to_state(class_finite_chain_state &state,const class_superblock &superblock);    /*!< Update the MPO stored at current position.*/
            extern void copy_superblock_env_to_state(class_finite_chain_state &state,const class_superblock &superblock);    /*!< Update the ENV stored at current position.*/
            extern int  insert_superblock_to_state(class_finite_chain_state &state, class_superblock &superblock);          /*!< Store current MPS and environments indexed by their respective positions on the chain. */
            extern int  move_center_point           (class_finite_chain_state & state , class_superblock & superblock);          /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
            extern void copy_state_to_superblock    (const class_finite_chain_state & state , class_superblock & superblock);    /*!< Update the MPS stored at current position.*/
        }


        namespace Ops {
            extern std::list<Eigen::Tensor<std::complex<double>,4>>
                        make_mpo_list                 (const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_L, const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_R);
            extern void apply_mpo                     (class_finite_chain_state &state,const Eigen::Tensor<std::complex<double>,4> mpo, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern void apply_mpo2                    (class_finite_chain_state &state,const Eigen::Tensor<std::complex<double>,4> mpo, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern void apply_mpo3                    (class_finite_chain_state &state,const Eigen::Tensor<std::complex<double>,4> mpo, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern void apply_mpos                    (class_finite_chain_state &state, const std::list<Eigen::Tensor<std::complex<double>,4>> &mpos, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern void normalize_chain               (class_finite_chain_state &state);
            extern void apply_energy_mpo_test         (class_finite_chain_state &state, class_superblock & superblock);
            extern void rebuild_environments          (class_finite_chain_state &state);
            extern void rebuild_superblock            (class_finite_chain_state &state, class_superblock & superblock);
            extern double overlap                     (const class_finite_chain_state &state1, const class_finite_chain_state &state2);
            extern double expectation_value           (const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>> &mpos, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern double exp_sq_value                (const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>> &mpos, const Eigen::Tensor<std::complex<double>,4> Ledge, const Eigen::Tensor<std::complex<double>,4> Redge);
            extern class_finite_chain_state
                        get_parity_projected_state    (const class_finite_chain_state &state, const Eigen::MatrixXcd paulimatrix, const int sign);
        }


        namespace Measure{

//            extern void do_all_measurements                           (class_finite_chain_state & state);
            extern int length                                      (const class_finite_chain_state & state);
            extern std::vector<int> bond_dimensions                   (const class_finite_chain_state & state);
            extern double norm                                        (const class_finite_chain_state & state);
            extern double energy_mpo                                  (class_finite_chain_state & state);
            extern double energy_per_site_mpo                         (class_finite_chain_state & state);
            extern double energy_variance_mpo                         (class_finite_chain_state & state);
            extern double energy_variance_per_site_mpo                (class_finite_chain_state & state);
            extern double midchain_entanglement_entropy               (const class_finite_chain_state & state);
            extern double spin_component(const class_finite_chain_state &state, const Eigen::Matrix2cd paulimatrix);
            extern Eigen::Tensor<std::complex<double>,1> mps_wavefn   (const class_finite_chain_state & state);
            extern std::vector<double> entanglement_entropies         (const class_finite_chain_state & state);
            extern std::vector<double> parities                       (class_finite_chain_state & state);

        }


        namespace Print {
            extern void print_full_state    (const class_finite_chain_state & state);
            extern void print_state         (const class_finite_chain_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_state_compact (const class_finite_chain_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_hamiltonians  (const class_finite_chain_state & state);
        }
        namespace Hdf5{
            extern void write_all_state                    (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_bond_matrices                (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_mps                    (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_mpo                    (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_env                    (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_env2                   (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_full_mps                     (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_full_mpo                     (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_hamiltonian_params           (const class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_entanglement                 (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_all_measurements             (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_all_parity_projections       (const class_finite_chain_state & state, const class_superblock &superblock, class_hdf5_file & hdf5, std::string sim_name);
            extern double write_parity_projected_analysis    (const class_finite_chain_state & state, const class_superblock &superblock, class_hdf5_file & hdf5, std::string sim_name,  std::string projection_name, const Eigen::MatrixXcd paulimatrix, const int sign);
            extern void load_from_hdf5                     (class_finite_chain_state & state, class_superblock &superblock, class_simulation_state &sim_state, class_hdf5_file &hdf5, std::string sim_name);
            extern void load_state_from_hdf5               (class_finite_chain_state & state, class_hdf5_file &hdf5, std::string sim_name);
            extern void load_sim_state_from_hdf5           (class_simulation_state &sim_state, class_hdf5_file &hdf5, std::string sim_name);
//            extern void load_sim_state_from_hdf5           (class_simulation_state &sim_state, class_hdf5_file &hdf5, std::string sim_name);

        }

        namespace Debug {
            extern void check_integrity          (const class_finite_chain_state &state, const class_superblock & superblock, class_simulation_state &sim_state);
            extern void check_integrity_of_sim   (const class_finite_chain_state &state, const class_superblock &superblock, class_simulation_state &sim_state);
            extern void check_integrity_of_mps   (const class_finite_chain_state &state);
            extern void print_parity_properties  (const class_finite_chain_state &state);

        }

    }




    namespace Infinite
    /*!
     * Functions for infinite MPS algorithms like iDMRG and iTEBD.
     * These functions do not require the class_finite_chain_state class, only the
     * class_superblock containing the local 2-site MPS + ENV representation.
     */
    {

        namespace Measure{

            namespace Results {
                extern int length;
                extern int bond_dimension;
                extern double norm;
                extern double truncation_error;
                extern double energy_mpo;
                extern double energy_ham;
                extern double energy_mom;
                extern double energy_per_site;
                extern double energy_variance;
                extern double energy_variance_per_site;
                extern double mps_wavefn;
                extern double midchain_entanglement_entropy;
                extern bool   superblock_measured;
            }
            extern void   do_all_measurements             (const class_superblock & superblock);
            extern int    length                          (const class_superblock & superblock);
            extern double norm                            (const class_superblock & superblock);
            extern double truncation_error                (const class_superblock & superblock);
            extern double energy_mpo                      (const class_superblock & superblock);
            extern double energy_ham                      (const class_superblock & superblock);
            extern double energy_mom                      (const class_superblock & superblock);
            extern double energy_variance_mpo             (const class_superblock & superblock);
            extern double energy_variance_ham             (const class_superblock & superblock);
            extern double energy_variance_mom             (const class_superblock & superblock);
            extern double midchain_entanglement_entropy   (const class_superblock & superblock);
        }

        namespace Print {
            extern void print_state         (const class_superblock & superblock);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_state_compact (const class_superblock & superblock);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_hamiltonians  (const class_superblock & superblock);
        }

        namespace Hdf5{
            extern void write_all_superblock(class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_2site_mps                    (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_2site_mpo                    (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_2site_env                    (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_2site_env2                   (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_hamiltonian_params           (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_all_measurements             (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void load_from_hdf5                     (class_superblock &superblock, class_simulation_state &sim_state, class_hdf5_file &hdf5, std::string sim_name);
            extern void load_sim_state_from_hdf5           (class_simulation_state &sim_state, class_hdf5_file &hdf5, std::string sim_name);
            extern void load_superblock_from_hdf5          (class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);

        }
    }





    namespace Common{

        namespace Info{

        }

        namespace Hdf5 {
            extern void write_simulation_state   (const class_simulation_state & sim_state,class_hdf5_file &hdf5,std::string sim_name);
        }


        namespace Measure {

            extern void   set_not_measured                (class_superblock & superblock);
            extern int    length                          (class_superblock & superblock);
            extern int    bond_dimension                  (class_superblock & superblock);
            extern double truncation_error                (class_superblock & superblock);
            extern double norm                            (const class_superblock & superblock);
            extern double energy_mpo                      (class_superblock & superblock);
            extern double energy_mpo                      (const class_superblock & superblock, const Eigen::Tensor<std::complex<double>,4> &theta);
            extern double energy_per_site_mpo             (class_superblock & superblock);
            extern double energy_per_site_ham             (class_superblock & superblock);
            extern double energy_per_site_mom             (class_superblock & superblock);
            extern double energy_variance_mpo             (class_superblock & superblock);
            extern double energy_variance_mpo             (const class_superblock & superblock, const Eigen::Tensor<std::complex<double>,4> &theta, double & energy_mpo);
            extern double energy_variance_mpo             (class_superblock & superblock, double &energy_mpo);
            extern double energy_variance_per_site_mpo    (class_superblock & superblock);
            extern double energy_variance_per_site_mpo    (class_superblock & superblock, double &energy_mpo);
            extern double energy_variance_per_site_ham    (class_superblock & superblock);
            extern double energy_variance_per_site_ham    (class_superblock & superblock, double &energy_per_site_ham);
            extern double energy_variance_per_site_mom    (class_superblock & superblock);
            extern double current_entanglement_entropy    (class_superblock & superblock);

        }



        namespace Views {
            extern Eigen::Tensor<std::complex<double>,4> theta, theta_evn_normalized, theta_odd_normalized;
            extern Eigen::Tensor<std::complex<double>,4> theta_sw ;
            extern Eigen::Tensor<std::complex<double>,3> LBGA, LAGB;
            extern Eigen::Tensor<std::complex<double>,2> l_evn, r_evn;
            extern Eigen::Tensor<std::complex<double>,2> l_odd, r_odd;
            extern Eigen::Tensor<std::complex<double>,4> transfer_matrix_LBGA;
            extern Eigen::Tensor<std::complex<double>,4> transfer_matrix_LAGB;
            extern Eigen::Tensor<std::complex<double>,4> transfer_matrix_evn;
            extern Eigen::Tensor<std::complex<double>,4> transfer_matrix_odd;
            extern void compute_mps_components(const class_superblock &superblock);

            extern Eigen::Tensor<std::complex<double>,4> get_theta                       (const class_finite_chain_state & state, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/

            extern Eigen::Tensor<std::complex<double>,4> get_theta                       (const class_superblock & superblock, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<std::complex<double>,4> get_theta_swapped               (const class_superblock & superblock, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<std::complex<double>,4> get_theta_evn                   (const class_superblock & superblock, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<std::complex<double>,4> get_theta_odd                   (const class_superblock & superblock, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_zero        (const class_superblock & superblock);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_LBGA        (const class_superblock & superblock, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_GALC        (const class_superblock & superblock, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_GBLB        (const class_superblock & superblock, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_LCGB        (const class_superblock & superblock, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_theta_evn   (const class_superblock & superblock, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_theta_odd   (const class_superblock & superblock, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_AB          (const class_superblock & superblock, int p);

            extern Eigen::Tensor<std::complex<double>,4> get_theta                       (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<std::complex<double>,4> get_theta_swapped               (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<std::complex<double>,4> get_theta_evn                   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<std::complex<double>,4> get_theta_odd                   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_zero        (const class_mps_2site  &MPS);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_LBGA        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_GALC        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_GBLB        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_LCGB        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_theta_evn   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_theta_odd   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<std::complex<double>,4> get_transfer_matrix_AB          (const class_mps_2site  &MPS, int p);


        }
    }

}








#endif //DMRG_NMSPC_FINITE_CHAIN_TOOLS_H
