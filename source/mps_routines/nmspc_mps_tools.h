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


namespace MPS_Tools{


    namespace Finite
    /*!
     * Functions for finite MPS algorithms like fDMRG and xDMRG
     * These functions require the class_finite_chain_state containing the
     * MPS and environments for all sites of the lattice.
     */
    {
        namespace Chain {
            extern void copy_superblock_to_chain    (class_finite_chain_state & state , const class_superblock & superblock);    /*!< Update the MPS, MPO and ENV stored at current position.*/
            extern void copy_superblock_mps_to_chain(class_finite_chain_state & state , const class_superblock & superblock);    /*!< Update the MPS stored at current position.*/
            extern void copy_superblock_mpo_to_chain(class_finite_chain_state & state , const class_superblock & superblock);    /*!< Update the MPO stored at current position.*/
            extern void copy_superblock_env_to_chain(class_finite_chain_state & state , const class_superblock & superblock);    /*!< Update the ENV stored at current position.*/
            extern int  insert_superblock_to_chain  (class_finite_chain_state & state , class_superblock & superblock);          /*!< Store current MPS and environments indexed by their respective positions on the chain. */
            extern int  move_center_point           (class_finite_chain_state & state , class_superblock & superblock);          /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
            extern void copy_chain_to_superblock    (const class_finite_chain_state & state , class_superblock & superblock);    /*!< Update the MPS stored at current position.*/
        }


        namespace Ops {
            extern void apply_mpo               (class_finite_chain_state &state,const Eigen::Tensor<std::complex<double>,4> mpo, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern void normalize_chain         (class_finite_chain_state &state);
            extern void rebuild_environments    (class_finite_chain_state &state);
            extern void rebuild_superblock      (class_finite_chain_state &state, class_superblock & superblock);
            extern void set_parity_projected_mps(class_finite_chain_state &state, const Eigen::MatrixXcd paulimatrix);
        }


        namespace Measure{
            namespace Results {
                extern size_t length;
                extern size_t bond_dimension;
                extern double norm;
                extern double truncation_error;
                extern double energy_mpo;
                extern double energy_per_site_mpo;
                extern double energy_variance_mpo;
                extern double midchain_entanglement_entropy;
                extern std::vector<double> entanglement_entropy;
                extern double parity_sx;
                extern double parity_sy;
                extern double parity_sz;
                extern Eigen::Tensor<std::complex<double>,1> mps_wavefn;
                extern bool state_measured;

            }


            extern void do_all_measurements                           (const class_finite_chain_state & state);
            extern size_t length                                      (const class_finite_chain_state & state);
            extern double norm                                        (const class_finite_chain_state & state);
            extern size_t bond_dimension                              (const class_finite_chain_state & state);
            extern double energy_mpo                                  (const class_finite_chain_state & state);
            extern double energy_variance_mpo                         (const class_finite_chain_state & state);
            extern double energy_variance_mpo                         (const class_finite_chain_state & state, double energy);
            extern double midchain_entanglement_entropy               (const class_finite_chain_state & state);
            extern std::vector<double>
            entanglement_entropy                                      (const class_finite_chain_state & state);
            extern double parity                                      (const class_finite_chain_state & state,const Eigen::Matrix2cd  paulimatrix);
            extern Eigen::Tensor<std::complex<double>,1> mps_wavefn   (const class_finite_chain_state & state);

        }


        namespace Print {
            extern void print_state         (const class_finite_chain_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_state_compact (const class_finite_chain_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_hamiltonians  (const class_finite_chain_state & state);
        }
        namespace Hdf5{
            extern void write_all                          (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_bond_matrices                (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_mps                    (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_mpo                    (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_env                    (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_2site_env2                   (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_full_mps                     (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_full_mpo                     (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_hamiltonian_params           (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
            extern void write_entanglement                 (class_finite_chain_state & state, class_hdf5_file & hdf5, std::string sim_name);
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
                extern size_t length;
                extern size_t bond_dimension;
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
            extern size_t length                          (const class_superblock & superblock);
            extern double norm                            (const class_superblock & superblock);
            extern size_t truncation_error                (const class_superblock & superblock);
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
            extern void write_all                          (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_bond_matrix                  (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_mps                          (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_mpo                          (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_env                          (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_env2                         (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_hamiltonian_params           (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
            extern void write_entanglement                 (const class_superblock &superblock, class_hdf5_file &hdf5, std::string sim_name);
        }
    }





    namespace Common{

        namespace Info{

        }




        namespace Measure {
            namespace Results {
                extern size_t length;
                extern size_t bond_dimension;
                extern double norm;
                extern double truncation_error;
                extern double energy_mpo, energy_per_site_mpo;
                extern double energy_ham;
                extern double energy_mom;
                extern double energy_variance_mpo, energy_variance_per_site_mpo;
                extern double energy_variance_ham;
                extern double energy_variance_mom;
                extern double midchain_entanglement_entropy;
                extern bool   superblock_measured;
            }

            extern void   set_not_measured();
            extern void   do_all_measurements             (const class_superblock & superblock);

            extern size_t length                          (const class_superblock & superblock);
            extern size_t bond_dimension                  (const class_superblock & superblock);
            extern double norm                            (const class_superblock & superblock);
            extern double truncation_error                (const class_superblock & superblock);

            extern double energy_mpo                      (const class_superblock & superblock);
            extern double energy_ham                      (const class_superblock & superblock);
            extern double energy_mom                      (const class_superblock & superblock);
            extern double energy_variance_mpo             (const class_superblock & superblock);
            extern double energy_variance_mpo             (const class_superblock & superblock, double & energy_mpo);
            extern double energy_variance_ham             (const class_superblock & superblock);
            extern double energy_variance_ham             (const class_superblock & superblock, double & energy_ham);
            extern double energy_variance_mom             (const class_superblock & superblock);
            extern double energy_variance_mom             (const class_superblock & superblock, double & energy_mom);
            extern double midchain_entanglement_entropy   (const class_superblock & superblock);

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
