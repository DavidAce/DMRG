//
// Created by david on 7/30/17.
//
#include <sim_parameters/nmspc_sim_settings.h>
#include <algorithms/class_algorithm_launcher.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_state.h>
#include <io/class_hdf5_file.h>
#include <io/class_hdf5_table_buffer2.h>
#include <algorithms/class_iDMRG.h>
#include <algorithms/class_fDMRG.h>
#include <algorithms/class_xDMRG.h>
#include <algorithms/class_iTEBD.h>
#include <spdlog/spdlog.h>


namespace s = settings;
using namespace std;
using namespace Textra;

class_algorithm_launcher::class_algorithm_launcher(std::shared_ptr<class_hdf5_file> hdf5_):hdf5(std::move(hdf5_))
{
    hdf5_path = hdf5->get_file_path();

}
class_algorithm_launcher::class_algorithm_launcher()
{
    hdf5 = std::make_shared<class_hdf5_file>(
            settings::hdf5::output_filename,
            settings::hdf5::output_folder,
            settings::hdf5::overwrite_file_if_found,
            settings::hdf5::resume_from_file,
            settings::hdf5::create_dir_if_not_found
            );
    hdf5_path = hdf5->get_file_path();
}


void class_algorithm_launcher::run_algorithms(){
    run_iDMRG();
    run_fDMRG();
    run_xDMRG();
    run_iTEBD();
    spdlog::info("All simulations finished");
    bool OK = true;
    hdf5->write_dataset(OK, "/common/fileOK");
}


void class_algorithm_launcher::run_iDMRG(){
    if(settings::idmrg::on){
        class_iDMRG iDMRG(hdf5);
        iDMRG.run();
    }
}


void class_algorithm_launcher::run_fDMRG(){
    if(settings::fdmrg::on){
        class_fDMRG fDMRG(hdf5);
        fDMRG.run();
    }
}

void class_algorithm_launcher::run_xDMRG(){
    if(settings::xdmrg::on){
        class_xDMRG xDMRG(hdf5);
        xDMRG.run();
    }
}

void class_algorithm_launcher::run_iTEBD(){
    if(settings::itebd::on){
        class_iTEBD iTEBD(hdf5);
        iTEBD.run();
    }
}

