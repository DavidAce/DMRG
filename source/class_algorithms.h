//
// Created by david on 7/30/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
#define FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H


#include <class_superblock.h>
#include <class_storage.h>
#include <class_tic_toc.h>
#include <class_hdf5.h>
#include <n_settings.h>
#include <iostream>
#include <iomanip>
#include <n_data_containers.h>

class class_algorithms  {
private:
    class_hdf5 hdf5;

public:
    class_algorithms():
    hdf5(class_hdf5(settings::hdf5::filename, settings::hdf5::path, true)){};

    void iDMRG(class_superblock &superblock, class_storage &storage);
    void fDMRG(class_superblock &superblock, class_storage &storage);
    void iTEBD(class_superblock &superblock);
    void FES(class_superblock &superblock);


private:
    void write_to_file_model(class_superblock &superblock){
        std::string group_name = "/Model/";
        hdf5.create_group(group_name);
        hdf5.write_to_file(superblock.H.asTensor4, group_name + "H/asTensor4");
        hdf5.write_to_file(superblock.H.W, group_name + "H/W");
        hdf5.write_to_file(superblock.H.asTimeEvolution, group_name + "H/asTimeEvolution");
        hdf5.write_to_file(superblock.H.asMatrix, group_name + "H/asMatrix");
    }

    void write_to_file_DMRG(class_superblock &superblock, int length) {
        std::ostringstream ostr;
        ostr << std::setfill('0') << std::setw(4) << std::abs(length);
        std::string group_name = "/DMRG/L_" + ostr.str() + "/";
        hdf5.create_group(group_name, "Length", length);
        hdf5.write_to_file(superblock.MPS.GA, group_name + "MPS/GA");
        hdf5.write_to_file(superblock.MPS.GB, group_name + "MPS/GB");
        hdf5.write_to_file(superblock.MPS.LA, group_name + "MPS/LA");
        hdf5.write_to_file(superblock.MPS.LB, group_name + "MPS/LB");
        hdf5.write_to_file(superblock.Lblock.block, group_name + "ENV/Lblock");
        hdf5.write_to_file(superblock.Rblock.block, group_name + "ENV/Rblock");
        hdf5.write_to_file(superblock.MPS.get_energy(superblock.H.asTensor4), group_name + "Energy");
        hdf5.write_to_file(superblock.MPS.get_entropy(), group_name + "Entropy");
    }

    void write_to_file_FES(tebd_data_container &data, long max_chi){
        std::ostringstream ostr;
        ostr << std::setfill('0') << std::setw(3) << std::abs(max_chi);
        std::string group_name = "/FES/chi_" + ostr.str() + "/";
        hdf5.create_group(group_name, "Max Chi", max_chi);
        hdf5.write_to_file(max_chi          , group_name + "max_chi");
        hdf5.write_to_file(data.chi         , group_name + "chi");
        hdf5.write_to_file(data.time_step   , group_name + "delta_t");
        hdf5.write_to_file(data.wall_time   , group_name + "wall_time");
        hdf5.write_to_file(data.phys_time   , group_name + "phys_time");
        hdf5.write_to_file(data.energy      , group_name + "energy");
        hdf5.write_to_file(data.entropy     , group_name + "entropy");
        hdf5.write_to_file(data.trunc_error , group_name + "trunc_error");
    }




};


#endif //FINITE_DMRG_EIGEN_CLASS_ALGORITHMS_H
