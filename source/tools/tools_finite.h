//
// Created by david on 2019-10-13.
//

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <complex>

class class_finite_state;
class class_model_base;
class class_simulation_status;


namespace tools{
    namespace finite
        /*!
         * Functions for finite MPS algorithms like fDMRG and xDMRG
         * These functions require the class_finite_chain_state containing the
         * MPS and environments for all sites of the lattice.
         */
    {
        namespace mps {
            extern void initialize                          (class_finite_state & state, size_t length);
            extern void randomize                           (class_finite_state &state, const std::string &parity_sector = "random", int seed_state = -1, bool use_pauli_eigenstates = false, bool enumeration =  false);
            extern void normalize                           (class_finite_state & state, const std::vector<size_t> & bond_dimensions = {});
            extern void rebuild_environments                (class_finite_state & state);
            extern int  move_center_point                   (class_finite_state & state);          /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
            extern void project_to_closest_parity_sector    (class_finite_state & state, std::string paulistring,bool keep_bond_dimensions = false);

            namespace internals{
                inline bool seed_state_unused = true;
                extern void set_product_state_in_parity_sector_from_bitset(class_finite_state & state, const std::string &parity_sector, const int seed_state);
                extern void set_product_state_in_parity_sector_randomly(class_finite_state & state, const std::string &parity_sector);
                extern void set_product_state_randomly(class_finite_state & state,const std::string &parity_sector,bool use_pauli_eigenstates);
            }

        }

        namespace mpo {
            extern void initialize                 (class_finite_state & state, size_t length, std::string model_type);
            extern void randomize                  (class_finite_state & state, int seed_state = -1);
            extern void reduce_mpo_energy          (class_finite_state & state);
        }

        namespace ops {
            extern std::list<Eigen::Tensor<Scalar,4>>
            make_mpo_list                 (const std::list<std::unique_ptr<class_model_base>> & mpos_L, const std::list<std::unique_ptr<class_model_base>> & mpos_R);
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
//            enum class OptMode  {OVERLAP, VARIANCE};
//            enum class OptSpace {SUBSPACE,DIRECT};
//            enum class OptType  {REAL, CPLX};


            extern Eigen::Tensor<Scalar,3> find_excited_state(const class_finite_state & state, const class_simulation_status & sim_status, OptMode optMode, OptSpace optSpace, OptType optType);
            extern Eigen::Tensor<Scalar,4> find_ground_state (const class_finite_state & state, std::string ritz = "SR");
            extern void truncate_theta(Eigen::Tensor<Scalar,3> & theta, class_finite_state & state, long chi_, double SVDThreshold);
            extern void truncate_left (Eigen::Tensor<Scalar,3> & theta, class_finite_state & state, long chi_, double SVDThreshold);
            extern void truncate_right(Eigen::Tensor<Scalar,3> & theta, class_finite_state & state, long chi_, double SVDThreshold);
            extern void truncate_theta(Eigen::Tensor<Scalar,4> & theta, class_finite_state & state, long chi_, double SVDThreshold);
        }

        namespace multisite{
            extern Eigen::DSizes<long,3> get_dimensions  (const class_finite_state &state, const std::list<size_t> &list_of_sites);
            extern size_t                get_problem_size(const class_finite_state &state, const std::list<size_t> &list_of_sites);
            extern std::list<size_t>     generate_site_list(class_finite_state &state, const size_t threshold, const size_t max_sites);
        }





        namespace print {
            extern void print_full_state    (const class_finite_state & state);
            extern void print_state         (const class_finite_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_state_compact (const class_finite_state & state);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
            extern void print_hamiltonians  (const class_finite_state & state);
        }

        namespace io{
            namespace internals{
                inline bool make_extendable_dataset(const std::string & prefix_path);
            }
            extern void write_all_state                              (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_bond_matrices                          (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_bond_matrix                            (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_full_mps                               (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_full_mpo                               (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_model                                  (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_entanglement                           (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_all_measurements                       (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path);
            extern void write_projection_to_closest_parity_sector    (const class_finite_state & state, h5pp::File & h5ppFile, const std::string & prefix_path, std::string parity_sector,bool keep_bond_dimensions);
            extern void load_from_hdf5                               (const h5pp::File & h5ppFile, class_finite_state & state    , class_simulation_status & sim_status, const std::string & prefix_path);
            extern class_finite_state load_state_from_hdf5           (const h5pp::File & h5ppFile, const std::string & prefix_path);
        }


        namespace debug {
            extern void check_integrity             (const class_finite_state & state);
            extern void check_integrity_of_mps      (const class_finite_state & state);
            extern void check_integrity_of_mpo      (const class_finite_state & state);
            extern void check_normalization_routine (const class_finite_state & state);
            extern void print_parity_properties     (const class_finite_state & state);

        }

    }
}
