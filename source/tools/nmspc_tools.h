//
// Created by david on 2019-01-29.
//

#ifndef tools_H
#define tools_H
#include <memory>
#include <string>
#include <general/nmspc_tensor_extra.h>
#include <io/nmspc_logger.h>
#include <general/class_tic_toc.h>

class class_infinite_state;
class class_mps_2site;
class class_finite_state;
class class_model_base;
class class_simulation_status;



namespace h5pp{
    class File;
}




namespace tools{
    inline std::shared_ptr<spdlog::logger> log;
    namespace finite
    /*!
     * Functions for finite MPS algorithms like fDMRG and xDMRG
     * These functions require the class_finite_chain_state containing the
     * MPS and environments for all sites of the lattice.
     */
    {
        using Scalar = std::complex<double>;
        namespace mps {
            extern void initialize                          (class_finite_state & state, size_t length);
            extern void randomize                           (class_finite_state & state, const std::string &parity_sector = "none", int seed_state = -1);
            extern void normalize                           (class_finite_state & state, bool keep_bond_dimensions = false);
            extern void rebuild_environments                (class_finite_state & state);
            extern int  move_center_point                   (class_finite_state & state);          /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
            extern void project_to_closest_parity_sector    (class_finite_state & state, std::string paulistring,bool keep_bond_dimensions = false);

            namespace internals{
                inline bool seed_state_unused = true;
                extern void set_product_state_in_parity_sector_from_bitset(class_finite_state & state, const std::string &parity_sector, const int seed_state);
                extern void set_product_state_in_parity_sector_randomly(class_finite_state & state, const std::string &parity_sector, const int seed_state);
                extern void set_product_state_randomly(class_finite_state & state,const std::string &parity_sector,  const int seed_state);
            }

        }

        namespace mpo {
            extern void initialize                 (class_finite_state & state, size_t length, std::string model_type);
            extern void randomize                  (class_finite_state & state);
        }

        namespace ops {
            extern std::list<Eigen::Tensor<Scalar,4>>
                        make_mpo_list                 (const std::list<std::shared_ptr<class_model_base>> & mpos_L, const std::list<std::shared_ptr<class_model_base>> & mpos_R);
            extern void apply_mpo                     (class_finite_state & state, const Eigen::Tensor<Scalar,4> & mpo, const Eigen::Tensor<Scalar,3> &Ledge, const Eigen::Tensor<Scalar,3> & Redge);
            extern void apply_mpos                    (class_finite_state & state, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
            extern class_finite_state
            get_projection_to_parity_sector           (const class_finite_state & state, const Eigen::MatrixXcd & paulimatrix, int sign,bool keep_bond_dimensions = false);
            extern class_finite_state
            get_projection_to_closest_parity_sector   (const class_finite_state & state, const Eigen::MatrixXcd & paulimatrix,bool keep_bond_dimensions = false);
            extern class_finite_state
            get_projection_to_closest_parity_sector   (const class_finite_state & state, std::string parity_sector, bool keep_bond_dimensions = false);
            extern double overlap                     (const class_finite_state & state1, const class_finite_state & state2);
            extern double expectation_value           (const class_finite_state & state1, const class_finite_state & state2,const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
            extern double exp_sq_value                (const class_finite_state & state1, const class_finite_state & state2,const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,4> & Ledge, const Eigen::Tensor<Scalar,4> & Redge);
        }


        namespace opt{
            enum class OptMode  {OVERLAP, VARIANCE};
            enum class OptSpace {SUBSPACE,DIRECT};
            enum class OptType  {REAL, CPLX};
            extern Eigen::Tensor<Scalar,3> find_excited_state(const class_finite_state & state, const class_simulation_status & sim_status, OptMode optMode, OptSpace optSpace, OptType optType);
            extern Eigen::Tensor<Scalar,4> find_ground_state (const class_finite_state & state, std::string ritz = "SR");
            extern void truncate_theta(Eigen::Tensor<Scalar,3> & theta, class_finite_state & state, long chi_, double SVDThreshold);
            extern void truncate_left (Eigen::Tensor<Scalar,3> & theta, class_finite_state & state, long chi_, double SVDThreshold);
            extern void truncate_right(Eigen::Tensor<Scalar,3> & theta, class_finite_state & state, long chi_, double SVDThreshold);
            extern void truncate_theta(Eigen::Tensor<Scalar,4> & theta, class_finite_state & state, long chi_, double SVDThreshold);
        }

        namespace multisite{
            extern std::list<size_t>  generate_site_list(class_finite_state &state, long threshold);
        }


        namespace measure{

//            extern void do_all_measurements                           (class_finite_state & state);
            extern int length                                         (const class_finite_state & state);
            extern size_t bond_dimension_current                      (const class_finite_state & state);
            extern size_t bond_dimension_midchain                     (const class_finite_state & state);
            extern std::vector<size_t> bond_dimensions                (const class_finite_state & state);
            extern double norm                                        (const class_finite_state & state);
            extern double energy                                      (const class_finite_state & state);
            extern double energy_per_site                             (const class_finite_state & state);
            extern double energy_variance                             (const class_finite_state & state);
            extern double energy_variance_per_site                    (const class_finite_state & state);
            extern double spin_component                              (const class_finite_state & state, Eigen::Matrix2cd paulimatrix);
            extern Eigen::Tensor<Scalar,1> mps_wavefn                 (const class_finite_state & state);
            extern double entanglement_entropy_current                (const class_finite_state & state);
            extern double entanglement_entropy_midchain               (const class_finite_state & state);
            extern std::vector<double> entanglement_entropies         (const class_finite_state & state);
            extern std::vector<double> spin_components                (const class_finite_state & state);

            namespace accurate{
                extern double energy                                    (const class_finite_state & state);
                extern double energy_per_site                           (const class_finite_state & state);
                extern double energy_variance                           (const class_finite_state & state);
                extern double energy_variance_per_site                  (const class_finite_state & state);
            }
            namespace multisite{
                extern double energy                                  (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
                extern double energy_per_site                         (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
                extern double energy_variance                         (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
                extern double energy_variance_per_site                (const class_finite_state & state, const Eigen::Tensor<Scalar,3> & multitheta);
            }

        }


        namespace print {
            extern void print_full_state    (const class_finite_state & state);
            extern void print_state         (const class_finite_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_state_compact (const class_finite_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_hamiltonians  (const class_finite_state & state);
        }

        namespace io{
            extern void write_all_state                              (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_bond_matrices                          (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_bond_matrix                            (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_full_mps                               (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_full_mpo                               (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_hamiltonian_params                     (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_entanglement                           (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_all_measurements                       (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_projection_to_closest_parity_sector    (const class_finite_state & state, h5pp::File & h5ppFile, std::string sim_name, std::string parity_sector,bool keep_bond_dimensions);
            extern void load_from_hdf5                               (const h5pp::File & h5ppFile, class_finite_state & state    , class_simulation_status & sim_status, std::string sim_name);
            extern class_finite_state load_state_from_hdf5           (const h5pp::File & h5ppFile, std::string sim_name);
        }


        namespace debug {
            extern void check_integrity             (const class_finite_state & state);
            extern void check_integrity_of_mps      (const class_finite_state & state);
            extern void check_integrity_of_mpo      (const class_finite_state & state);
            extern void check_normalization_routine (const class_finite_state & state);
            extern void print_parity_properties     (const class_finite_state & state);

        }

    }




    namespace infinite
    /*!
     * Functions for infinite MPS algorithms like iDMRG and iTEBD.
     * These functions do not require the class_finite_chain_state class, only the
     * class_superblock containing the local 2-site MPS + ENV representation.
     */
    {
        using Scalar = std::complex<double>;

        namespace mps{
            extern class_infinite_state set_random_state(const class_infinite_state & state, std::string parity);
        }

        namespace opt{
            extern Eigen::Tensor<Scalar,4> find_ground_state(const class_infinite_state & state, std::string ritz = "SR");
            extern Eigen::Tensor<Scalar,4> time_evolve_theta(const class_infinite_state & state, const Eigen::Tensor<Scalar, 4> &U);
            extern void truncate_theta(Eigen::Tensor<Scalar,4> &theta, class_infinite_state & state, long chi_, double SVDThreshold);

        }

        namespace measure{
            extern int    length                          (const class_infinite_state & state);
            extern int    bond_dimension                  (const class_infinite_state & state);
            extern double truncation_error                (const class_infinite_state & state);
            extern double norm                            (const class_infinite_state & state);
            extern double energy_mpo                      (const class_infinite_state & state);
            extern double energy_mpo                      (const class_infinite_state & state, const Eigen::Tensor<Scalar,4> &theta);
            extern double energy_per_site_mpo             (const class_infinite_state & state);
            extern double energy_per_site_ham             (const class_infinite_state & state);
            extern double energy_per_site_mom             (const class_infinite_state & state);
            extern double energy_variance_mpo             (const class_infinite_state & state, const Eigen::Tensor<Scalar,4> &theta, double &energy_mpo);
            extern double energy_variance_mpo             (const class_infinite_state & state, const Eigen::Tensor<Scalar,4> &theta);
            extern double energy_variance_mpo             (const class_infinite_state & state);
            extern double energy_variance_per_site_mpo    (const class_infinite_state & state);
            extern double energy_variance_per_site_ham    (const class_infinite_state & state);
            extern double energy_variance_per_site_mom    (const class_infinite_state & state);
            extern double current_entanglement_entropy    (const class_infinite_state & state);
        }

        namespace print {
            extern void print_state         (const class_infinite_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_state_compact (const class_infinite_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_hamiltonians  (const class_infinite_state & state);
        }

        namespace io{
            extern void write_all_state(const class_infinite_state &state, h5pp::File &h5ppFile, std::string sim_name);
            extern void write_2site_mps                    (const class_infinite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_2site_mpo                    (const class_infinite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_2site_env                    (const class_infinite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_2site_env2                   (const class_infinite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_hamiltonian_params           (const class_infinite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void write_all_measurements             (const class_infinite_state & state, h5pp::File & h5ppFile, std::string sim_name);
            extern void load_from_hdf5                     (const h5pp::File & h5ppFile, class_infinite_state & state, class_simulation_status &sim_status, std::string sim_name);
            extern void load_superblock_from_hdf5          (const h5pp::File & h5ppFile, class_infinite_state & state, std::string sim_name);
            extern void load_sim_status_from_hdf5           (const h5pp::File & h5ppFile, class_simulation_status & sim_status, std::string sim_name);

        }



        namespace debug {
            extern void check_integrity             (const class_infinite_state & state, const class_simulation_status & sim_status);
            extern void check_integrity_of_sim      (const class_infinite_state & state, const class_simulation_status & sim_status);
            extern void check_integrity_of_mps      (const class_infinite_state & state);
            extern void check_normalization_routine (const class_infinite_state & state);

        }


    }





    namespace common{
        using Scalar = std::complex<double>;

        namespace io {
            extern void
            write_simulation_status(const class_simulation_status &sim_status, h5pp::File &h5ppFile,
                                    std::string sim_name);
            extern class_simulation_status
            load_sim_status_from_hdf5(const h5pp::File &h5ppFile, std::string sim_name);

        }


        namespace profile{
            inline class_tic_toc t_eig;
            inline class_tic_toc t_svd;
            inline class_tic_toc t_ene;
            inline class_tic_toc t_var;
            inline class_tic_toc t_ent;
            inline class_tic_toc t_hdf;
            inline class_tic_toc t_prj;
            inline class_tic_toc t_opt;
            inline class_tic_toc t_chk;

            inline class_tic_toc t_ene_mpo;
            inline class_tic_toc t_ene_ham;
            inline class_tic_toc t_ene_mom;
            inline class_tic_toc t_var_mpo;
            inline class_tic_toc t_var_ham;
            inline class_tic_toc t_var_mom;

            extern void print_profiling(class_tic_toc &t_parent);
            extern void init_profiling();
        }


        namespace views {
            extern Eigen::Tensor<Scalar,4> theta, theta_evn_normalized, theta_odd_normalized;
            extern Eigen::Tensor<Scalar,4> theta_sw ;
            extern Eigen::Tensor<Scalar,3> LBGA, LAGB;
            extern Eigen::Tensor<Scalar,2> l_evn, r_evn;
            extern Eigen::Tensor<Scalar,2> l_odd, r_odd;
            extern Eigen::Tensor<Scalar,4> transfer_matrix_LBGA;
            extern Eigen::Tensor<Scalar,4> transfer_matrix_LAGB;
            extern Eigen::Tensor<Scalar,4> transfer_matrix_evn;
            extern Eigen::Tensor<Scalar,4> transfer_matrix_odd;
            extern bool components_computed;
            extern void compute_mps_components(const class_infinite_state &state);

            extern Eigen::Tensor<Scalar,4> get_theta                       (const class_finite_state & state, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/

            extern Eigen::Tensor<Scalar,4> get_theta                       (const class_infinite_state & state, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<Scalar,4> get_theta_swapped               (const class_infinite_state & state, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<Scalar,4> get_theta_evn                   (const class_infinite_state & state, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<Scalar,4> get_theta_odd                   (const class_infinite_state & state, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_zero        (const class_infinite_state & state);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LBGA        (const class_infinite_state & state, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GALC        (const class_infinite_state & state, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GBLB        (const class_infinite_state & state, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LCGB        (const class_infinite_state & state, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_evn   (const class_infinite_state & state, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_odd   (const class_infinite_state & state, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_AB          (const class_infinite_state & state, int p);

            extern Eigen::Tensor<Scalar,4> get_theta                       (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<Scalar,4> get_theta_swapped               (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<Scalar,4> get_theta_evn                   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
            extern Eigen::Tensor<Scalar,4> get_theta_odd                   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_zero        (const class_mps_2site  &MPS);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LBGA        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GALC        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GBLB        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LCGB        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_evn   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_odd   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
            extern Eigen::Tensor<Scalar,4> get_transfer_matrix_AB          (const class_mps_2site  &MPS, int p);


        }
    }

}








#endif //DMRG_NMSPC_FINITE_CHAIN_TOOLS_H
