//
// Created by david on 7/30/17.
//

#include <class_algorithms.h>
#include <n_model.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
namespace s = settings;


void class_algorithms::iDMRG(class_superblock &superblock, class_storage &storage){
/*!
 * \fn void iDMRG(class_superblock &superblock, class_storage &S, int max_length)
 * \brief Infinite DMRG, grows the chain from 2 up to `max_idmrg::length` particles.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param max_length Maximum chain length after which the algorithm stops.
 */
    using namespace settings::profiling;
    class_tic_toc t_tot(on, precision, "iDMRG Total Time           ");
    class_tic_toc t_eig(on, precision, "iDMRG Eigenvalue solver    ");
    class_tic_toc t_svd(on, precision, "iDMRG SVD Truncation       ");
    class_tic_toc t_env(on, precision, "iDMRG Enlarge environment  ");
    class_tic_toc t_sto(on, precision, "iDMRG Store MPS            ");
    class_tic_toc t_mps(on, precision, "iDMRG Update MPS           ");



    storage.set_length(s::idmrg::max_length);
    t_tot.tic();
    int length = 0;
    while(length < s::idmrg::max_length){
                        length += 2;
                        superblock.update_bond_dimensions();
        t_eig.tic();    superblock.find_ground_state(s::precision::eigSteps, s::precision::eigThreshold);    t_eig.toc();
        t_svd.tic();    superblock.truncate         (s::idmrg::max_chi,      s::precision::SVDThreshold);    t_svd.toc();
        t_mps.tic();    superblock.update_MPS();                                                             t_mps.toc();
        t_sto.tic();    storage.store_insert(superblock);                                                    t_sto.toc();
        t_env.tic();    superblock.enlarge_environment();                                                    t_env.toc();
        superblock.print_picture(s::console::graphics);
                        superblock.print_state(s::console::verbosity);
                        write_to_file_DMRG(superblock,length);
                        superblock.swap_AB();
    }
    t_tot.toc();
    superblock.print_state(s::console::verbosity + 1);
    t_eig.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_env.print_time_w_percent();
    t_sto.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
    write_to_file_model(superblock);
}



void class_algorithms::fDMRG(class_superblock &superblock, class_storage &storage){
/*!
 * \fn void fDMRG(class_superblock &superblock, class_storage &S, int sweeps)
 * \brief Finite DMRG sweeps across the chain built during iDMRG.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param sweeps Maximum number of sweeps.
 */
    using namespace settings::profiling;
    class_tic_toc t_tot(on, precision, "fDMRG Total Time           ");
    class_tic_toc t_eig(on, precision, "fDMRG Eigenvalue solver    ");
    class_tic_toc t_svd(on, precision, "fDMRG SVD Truncation       ");
    class_tic_toc t_env(on, precision, "fDMRG Enlarge environment  ");
    class_tic_toc t_sto(on, precision, "fDMRG Store MPS            ");
    class_tic_toc t_mps(on, precision, "fDMRG Update MPS           ");
    int direction  = 1;
    int sweep = 0;
    t_tot.tic();
//    Eigen::ArrayXi chi_list = Eigen::ArrayXi::LinSpaced(s::fdmrg::max_sweeps,
//                                                        s::chi,
//                                                        s::increasing_chi?
//                                                        std::max(s::chi_max, s::chi+10)
//                                                        : s::chi);

    while(sweep < s::fdmrg::max_sweeps) {
//                        int chi = chi_list(sweep);
                        storage.load(superblock);
                        superblock.update_bond_dimensions();
        t_eig.tic();    superblock.find_ground_state(s::precision::eigSteps, s::precision::eigThreshold);         t_eig.toc();
        t_svd.tic();    superblock.truncate         (s::fdmrg::max_chi, s::precision::SVDThreshold);    t_svd.toc();
        t_mps.tic();    superblock.update_MPS();                                            t_mps.toc();
        t_sto.tic();    storage.overwrite_MPS(superblock);                                  t_sto.toc();
                        superblock.print_picture(s::console::graphics);
                        superblock.print_state(s::console::verbosity);
        t_env.tic();    superblock.enlarge_environment(direction);                          t_env.toc();
        t_sto.tic();    storage.move(superblock, direction, sweep);                         t_sto.toc();
    }

    t_tot.toc();
    superblock.print_state(s::console::verbosity + 1);
    t_eig.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_env.print_time_w_percent();
    t_sto.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
}


void class_algorithms::iTEBD(class_superblock &superblock){
/*!
 * \fn iTEBD(class_superblock &superblock, class_hdf5 &hdf5)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param max_iter Maximum number of iterations.
 */
    using namespace settings::profiling;
    class_tic_toc t_evo(on, precision, "iTEBD Time evolution       ");
    class_tic_toc t_svd(on, precision, "iTEBD SVD Truncation       ");
    class_tic_toc t_mps(on, precision, "iTEBD Update MPS           ");
    class_tic_toc t_upd(on, precision, "iTEBD Update Evolution     ");
    class_tic_toc t_tot(on, precision, "iTEBD Total Time           ");
    superblock.reset();
    t_tot.tic();

    double delta_t = s::itebd::delta_t;
    int chi = settings::itebd::max_chi;
    double old_energy = 0;
    double new_energy= superblock.MPS.get_energy(superblock.H.asTensor4)(0);
    for(auto step = 0; step < s::itebd::max_steps ; step++){
                        superblock.update_bond_dimensions();
        t_evo.tic();    superblock.time_evolve();                                   t_evo.toc();
        t_svd.tic();    superblock.truncate   (chi, s::precision::SVDThreshold);    t_svd.toc();
        t_mps.tic();    superblock.update_MPS();                                    t_mps.toc();
                        superblock.swap_AB();

        t_upd.tic();
        if (Math::mod(step,500) == 0) {
            superblock.print_state(1, "itebd:");
//            cout << "step: " << step << " delta = " <<setprecision(16)<< delta_t
//                 << " diff: " << std::fabs((new_energy-old_energy)/old_energy) <<endl;
        }
        if(Math::mod(step,1)==0 && delta_t > 1e-8){
            old_energy = new_energy;
            new_energy = superblock.MPS.get_energy(superblock.H.asTensor4)(0);
            if (std::fabs((new_energy-old_energy)/old_energy) < 1e-6){
                delta_t *=0.1;
                superblock.H.asTimeEvolution = Model::TimeEvolution_4th_order(2,delta_t);
            }
        }
        t_upd.toc();
    }
    t_tot.toc();
    superblock.print_state(s::console::verbosity+1);
    t_evo.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_upd.print_time_w_percent();
    t_tot.print_time();
}

void class_algorithms::FES(class_superblock &superblock){
/*!
 * \fn FES(class_superblock &superblock,class_storage& storage, c  class_hdf5 &hdf5)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param max_iter Maximum number of iterations.
 */
    using namespace settings::profiling;
    class_tic_toc t_evo(on, precision, "FES Time evolution       ");
    class_tic_toc t_svd(on, precision, "FES SVD Truncation       ");
    class_tic_toc t_mps(on, precision, "FES Update MPS           ");
    class_tic_toc t_upd(on, precision, "FES Update Evolution     ");
    class_tic_toc t_tot(on, precision, "FES Total Time           ");

    t_tot.tic();

    Eigen::ArrayXi chi_list = Eigen::ArrayXi::LinSpaced(s::fes::num_chi,
                                                        s::fes::min_chi,
                                                        s::fes::max_chi);
    tebd_data_container data;
    data.clear();

    for(auto chi = chi_list.data() ; chi != chi_list.data() + chi_list.size() ; ++chi) {
        superblock.reset();
        double delta_t = s::itebd::delta_t;
        double phys_time = 0;
        int step = 0;

        superblock.H.asTimeEvolution = Model::TimeEvolution_4th_order(2,delta_t);
        double old_energy = 0;
        double new_energy= superblock.MPS.get_energy(superblock.H.asTensor4)(0);
        double old_entropy = 0;
        double new_entropy= superblock.MPS.get_entropy()(0);
        while(delta_t > 1e-4){
//        for (auto step = 0; step < s::fes::max_steps; step++) {
            superblock.update_bond_dimensions();
            t_evo.tic();    superblock.time_evolve();                               t_evo.toc();
            t_svd.tic();    superblock.truncate   (*chi, s::precision::SVDThreshold);           t_svd.toc();
            t_mps.tic();    superblock.update_MPS();                                t_mps.toc();
                            superblock.swap_AB();

            data.push_back(
                    superblock.MPS.GB.dimension(2),
                    delta_t,
                    phys_time,
                    t_tot.get_measured_time(),
                    superblock.MPS.get_energy(superblock.H.asTensor4)(0),
                    superblock.MPS.get_entropy()(0),
                    superblock.truncation_error
            );

            t_upd.tic();

            if(Math::mod(step,1)==0){
//                old_energy = new_energy;
//                new_energy = superblock.MPS.get_energy(superblock.H.asTensor4)(0);
                old_entropy = new_entropy;
                new_entropy= superblock.MPS.get_entropy()(0);
//                if (std::fabs((new_energy-old_energy)/old_energy)/delta_t < 1e-8){
                if (std::fabs((new_entropy-old_entropy)/old_entropy) < 1e-8){

                    cout << "hej" << endl;
                    delta_t *=0.1;
                    superblock.H.asTimeEvolution = Model::TimeEvolution_4th_order(2,delta_t);
                }
            }
            if (Math::mod(step,1000) == 0) {
                superblock.print_state(2, "itebd:");
                cout << "step: " << step << " delta = " <<setprecision(16)<< delta_t
                     //                        << " diff: " << std::fabs((new_energy-old_energy)/new_energy)/delta_t
                     << " diff: " << std::fabs((new_entropy-old_entropy)/old_entropy)
                     << " chi: " << *chi << '\n';
            }
            t_upd.toc();


            phys_time += delta_t;
            step++;
        }
        superblock.print_state(s::console::verbosity+1);
        write_to_file_FES(data, *chi);
        data.clear();
    }

    t_tot.toc();
    t_evo.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_upd.print_time_w_percent();
    t_tot.print_time();
}


