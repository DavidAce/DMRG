//
// Created by david on 2019-01-29.
//

#ifndef MPS_TOOLS_H
#define MPS_TOOLS_H
#include <memory>
#include <string>
#include <general/nmspc_tensor_extra.h>
#include <io/nmspc_logger.h>
#include <general/class_tic_toc.h>

class class_superblock;
class class_mps_2site;
class class_finite_chain_state;
class class_hdf5_file;
class class_hamiltonian_base;
class class_simulation_state;



namespace h5pp{
    class File;
}




namespace MPS_Tools{
    inline std::shared_ptr<spdlog::logger> log;
    namespace Finite
    /*!
     * Functions for finite MPS algorithms like fDMRG and xDMRG
     * These functions require the class_finite_chain_state containing the
     * MPS and environments for all sites of the lattice.
     */
    {
        namespace Chain {
            extern void initialize_state(class_finite_chain_state &state, std::string model_type, std::string parity, const size_t length);
            extern void randomize_mpos  (class_finite_chain_state &state);

            extern void copy_superblock_to_state        (class_finite_chain_state &state, const class_superblock & superblock);    /*!< Update the MPS, MPO and ENV stored at current position.*/
            extern void copy_superblock_mps_to_state    (class_finite_chain_state &state, const class_superblock & superblock);    /*!< Update the MPS stored at current position.*/
            extern void copy_superblock_mpo_to_state    (class_finite_chain_state &state, const class_superblock & superblock);    /*!< Update the MPO stored at current position.*/
            extern void copy_superblock_env_to_state    (class_finite_chain_state &state, const class_superblock & superblock);    /*!< Update the ENV stored at current position.*/
            extern int  insert_superblock_to_state      (class_finite_chain_state &state, const class_superblock & superblock);          /*!< Store current MPS and environments indexed by their respective positions on the chain. */
            extern int  move_center_point               (class_finite_chain_state &state, class_superblock & superblock);          /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
            extern void copy_state_to_superblock        (const class_finite_chain_state & state , class_superblock & superblock);    /*!< Update the MPS stored at current position.*/
        }


        namespace Ops {
            extern std::list<Eigen::Tensor<std::complex<double>,4>>
                        make_mpo_list                 (const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_L, const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_R);
            extern void apply_mpo                     (class_finite_chain_state &state,const Eigen::Tensor<std::complex<double>,4> mpo, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern void apply_mpos                    (class_finite_chain_state &state, const std::list<Eigen::Tensor<std::complex<double>,4>> &mpos, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
            extern void normalize_chain               (class_finite_chain_state &state);

            extern double exp_sq_value                (const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>> &mpos, const Eigen::Tensor<std::complex<double>,4> Ledge, const Eigen::Tensor<std::complex<double>,4> Redge);
            extern class_finite_chain_state
            set_random_product_state                  (const class_finite_chain_state &state, const std::string parity);
            extern class_finite_chain_state
            get_parity_projected_state                (const class_finite_chain_state &state, const Eigen::MatrixXcd paulimatrix, const int sign);
            extern class_finite_chain_state
            get_closest_parity_state                  (const class_finite_chain_state &state, const Eigen::MatrixXcd paulimatrix);
            extern class_finite_chain_state
            get_closest_parity_state                  (const class_finite_chain_state &state, const std::string paulistring);
            extern void rebuild_environments          (class_finite_chain_state &state);
            extern void rebuild_superblock            (const class_finite_chain_state &state, class_superblock & superblock);
            extern double overlap                     (const class_finite_chain_state &state1, const class_finite_chain_state &state2);
            extern double expectation_value           (const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>> &mpos, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge);
        }


        namespace Opt{
            enum class OptMode  {OVERLAP, VARIANCE};
            enum class OptSpace {PARTIAL,FULL,DIRECT,GUIDED,CPPOPTLIB};
            extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> find_optimal_excited_state(const class_superblock & superblock, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace);

            namespace internals{
                extern std::tuple<Eigen::MatrixXd, Eigen::VectorXd>              find_subspace          (const class_superblock & superblock, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace);
                extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> direct_optimization    (const class_superblock & superblock, const class_simulation_state & sim_state);
                extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> guided_optimization    (const class_superblock & superblock, const class_simulation_state & sim_state);
                extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> cppoptlib_optimization (const class_superblock & superblock, const class_simulation_state & sim_state);
                extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> subspace_optimization  (const class_superblock & superblock, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace);
                extern std::vector<int> generate_size_list(const int shape);
                class base_functor;
                class subspace_functor;
                class direct_functor;
            }

        }

        namespace Measure{

//            extern void do_all_measurements                           (class_finite_chain_state & state);
            extern int length                                         (const class_finite_chain_state & state);
            extern std::vector<int> bond_dimensions                   (const class_finite_chain_state & state);
            extern double norm                                        (const class_finite_chain_state & state);
            extern double energy_mpo                                  (const class_finite_chain_state & state);
            extern double energy_per_site_mpo                         (const class_finite_chain_state & state);
            extern double energy_variance_mpo                         (const class_finite_chain_state & state);
            extern double energy_variance_per_site_mpo                (const class_finite_chain_state & state);
            extern double midchain_entanglement_entropy               (const class_finite_chain_state & state);
            extern double spin_component                              (const class_finite_chain_state &state, const Eigen::Matrix2cd paulimatrix);
            extern Eigen::Tensor<std::complex<double>,1> mps_wavefn   (const class_finite_chain_state & state);
            extern std::vector<double> entanglement_entropies         (const class_finite_chain_state & state);
            extern std::vector<double> spin_components                (const class_finite_chain_state & state);
        }


        namespace Print {
            extern void print_full_state    (const class_finite_chain_state & state);
            extern void print_state         (const class_finite_chain_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_state_compact (const class_finite_chain_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_hamiltonians  (const class_finite_chain_state & state);
        }

        namespace H5pp{
            extern void write_all_state                    (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_bond_matrices                (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_bond_matrix                  (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_full_mps                     (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_full_mpo                     (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_hamiltonian_params           (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_entanglement                 (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_all_measurements             (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_closest_parity_projection    (const class_finite_chain_state & state, h5pp::File & h5ppFile, std::string sim_name, std::string paulistring);
            extern void load_from_hdf5                     (const h5pp::File & h5ppFile, class_finite_chain_state & state    , class_superblock & superblock, class_simulation_state & sim_state, std::string sim_name);
            extern void load_state_from_hdf5               (const h5pp::File & h5ppFile, class_finite_chain_state & state    , std::string sim_name);
            extern void load_sim_state_from_hdf5           (const h5pp::File & h5ppFile, class_simulation_state   & sim_state, std::string sim_name);
//            extern void load_sim_state_from_hdf5           (class_simulation_state &sim_state, class_hdf5_file &hdf5, std::string sim_name);

        }

        namespace Debug {
            extern void check_integrity             (const class_finite_chain_state & state, const class_superblock & superblock, const class_simulation_state & sim_state);
            extern void check_integrity_of_sim      (const class_finite_chain_state & state, const class_superblock & superblock, const class_simulation_state & sim_state);
            extern void check_integrity_of_mps      (const class_finite_chain_state & state);
            extern void check_normalization_routine (const class_finite_chain_state & state);
            extern void print_parity_properties     (const class_finite_chain_state & state);

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

        namespace H5pp{
            extern void write_all_superblock               (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_2site_mps                    (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_2site_mpo                    (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_2site_env                    (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_2site_env2                   (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_hamiltonian_params           (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_all_measurements             (const class_superblock & superblock, h5pp::File & h5ppFile, std::string sim_name);
            extern void load_from_hdf5                     (const h5pp::File & h5ppFile, class_superblock & superblock, class_simulation_state &sim_state, std::string sim_name);
            extern void load_superblock_from_hdf5          (const h5pp::File & h5ppFile, class_superblock & superblock, std::string sim_name);
            extern void load_sim_state_from_hdf5           (const h5pp::File & h5ppFile, class_simulation_state & sim_state, std::string sim_name);

        }

    }





    namespace Common{

        namespace Info{

        }

        namespace Prof{
            namespace Obs{
                inline class_tic_toc t_eig;
                inline class_tic_toc t_ene_mpo;
                inline class_tic_toc t_ene_ham;
                inline class_tic_toc t_ene_mom;
                inline class_tic_toc t_var_mpo;
                inline class_tic_toc t_var_ham;
                inline class_tic_toc t_var_mom;
                inline class_tic_toc t_entropy;
                inline class_tic_toc t_temp1  ;
                inline class_tic_toc t_temp2  ;
                inline class_tic_toc t_temp3  ;
                inline class_tic_toc t_temp4  ;
                extern void print_profiling(class_tic_toc &t_parent);
            }
            extern void enable_profiling(int precision = 5);

        }

        namespace H5pp {
            extern void write_algorithm_state(const class_simulation_state &sim_state, h5pp::File &h5ppFile,
                                              std::string sim_name);
        }

        namespace Measure {
            extern int    length                          (const class_superblock & superblock);
            extern int    bond_dimension                  (const class_superblock & superblock);
            extern double truncation_error                (const class_superblock & superblock);
            extern double norm                            (const class_superblock & superblock);
            extern double energy_mpo                      (const class_superblock & superblock);
            extern double energy_mpo                      (const class_superblock & superblock, const Eigen::Tensor<std::complex<double>,4> &theta);
            extern double energy_per_site_mpo             (const class_superblock & superblock);
            extern double energy_per_site_ham             (const class_superblock & superblock);
            extern double energy_per_site_mom             (const class_superblock & superblock);
            extern double energy_variance_mpo             (const class_superblock & superblock, const Eigen::Tensor<std::complex<double>,4> &theta, double &energy_mpo);
            extern double energy_variance_mpo             (const class_superblock & superblock, const Eigen::Tensor<std::complex<double>,4> &theta);
            extern double energy_variance_mpo             (const class_superblock & superblock);
            extern double energy_variance_per_site_mpo    (const class_superblock & superblock);
            extern double energy_variance_per_site_mpo    (const class_superblock & superblock);
            extern double energy_variance_per_site_ham    (const class_superblock & superblock);
            extern double energy_variance_per_site_ham    (const class_superblock & superblock);
            extern double energy_variance_per_site_mom    (const class_superblock & superblock);
            extern double current_entanglement_entropy    (const class_superblock & superblock);
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
            extern bool components_computed;
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
