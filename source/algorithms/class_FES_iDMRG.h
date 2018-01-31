//
// Created by david on 2018-01-31.
//

#ifndef DMRG_CLASS_FES_DMRG_H
#define DMRG_CLASS_FES_DMRG_H



#include "class_base.h"
/*!
 *
 * \brief Class that runs the Finite-entanglement scaling (FES) mode in the DMRG algorithm.
 *
 * # Finite-entanglement scaling (FES) DMRG class
 * \param shared_ptr<class_hdf5_file> An hdf5 class object that handles the output file.
 * \param shared_ptr<class_hdf5_table_buffer> (optional) A buffer for table entries that goes into the output file
 * \param shared_ptr<class_superblock> (optional) A class that stores current MPS and environments at each iteration.
 * \param shared_ptr<class_measurement> (optional) A class that extracts, or measures, quantities from the superblock.
 */
class class_infinite_DMRG;
class class_FES_iDMRG : public class_base {
public:
    //Inherit the constructor of class_base
    using class_base::class_base;
    explicit class_FES_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    long chi_min    = settings::fes_idmrg::chi_min;
    long chi_max    = settings::fes_idmrg::chi_max;
    long chi_num    = settings::fes_idmrg::chi_num;
    int  max_length = settings::fes_idmrg::max_length;
    int  print_freq = settings::fes_idmrg::print_freq;

    int  iteration  = 0;

    void run() override;
    void run2();

    void print_status_full()   override;
    void print_status_update() override;

    void store_table_entry()   override;
    void print_profiling()     override;
};





#endif //DMRG_CLASS_FES_DMRG_H
