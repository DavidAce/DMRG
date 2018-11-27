//
// Created by david on 2018-02-09.
//

#ifndef DMRG_CLASS_EXITED_DMRG_H
#define DMRG_CLASS_EXITED_DMRG_H

#include "class_algorithm_base.h"
#include <unsupported/Eigen/CXX11/Tensor>
class class_table_finite_chain;
class class_table_dmrg;

/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_finite_chain_sweeper;
class class_xDMRG : public class_algorithm_base {
private:

public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_xdmrg;
    std::unique_ptr<class_hdf5_table<class_table_finite_chain>> table_xdmrg_chain;

    enum class xDMRG_Mode {FULL,PARTIAL};
    int    full_iters   = 5;
    int    max_length   ;
    int    max_sweeps   ;

    //Energy ranges
    double energy_min = 0;
    double energy_max = 0;
    double energy_target = 0;
    double energy_now = 0;

    void run()                                          override;
    void check_convergence_all()                        override;
    void initialize_constants()                         override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_table_entry_to_file()                    override;
    void store_chain_entry_to_file();
    void store_mps_to_file(bool force = false);
    void single_xDMRG_step();
    void initialize_chain();
    void reset_chain_mps_to_random_product_state(std::string parity = "none");
    void set_random_fields_in_chain_mpo();
    void find_energy_range();

    std::vector<int> generate_size_list(const int shape);

    template <typename Derived>
    std::vector<size_t> make_sorted_index(const Eigen::MatrixBase<Derived>& values)
    {
        std::vector<size_t> index(values.derived().size());
        std::iota(index.begin(), index.end(), 0);
        std::sort(index.begin(), index.end(), [&values](size_t a, size_t b) { return std::real(values[a]) > std::real(values[b]); } );
        return index;
    }

    void sort_and_filter_eigenstates(Eigen::VectorXcd &eigvals,
                                     Eigen::MatrixXcd &eigvecs,
                                     Eigen::VectorXd  &overlaps,
                                     int &nev,
                                     double overlap_cutoff = 1e-2);
    Eigen::Tensor<Scalar,4> find_state_with_greatest_overlap_full_diag (Eigen::Tensor<Scalar, 4> &theta);
    Eigen::Tensor<Scalar,4> find_state_with_greatest_overlap_part_diag3 (Eigen::Tensor<Scalar, 4> &theta, xDMRG_Mode mode = xDMRG_Mode::PARTIAL);
    Eigen::Tensor<Scalar,4> find_state_with_greatest_state_in_subspace (Eigen::Tensor<Scalar, 4> &theta);



};



#endif //DMRG_CLASS_EXITED_DMRG_H
