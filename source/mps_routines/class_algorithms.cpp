//
// Created by david on 7/30/17.
//
//#define EIGEN_USE_MKL_ALL
#include <mps_routines/class_algorithms.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_sweep_storage.h>
#include <mps_routines/class_observables.h>
#include <IO/class_multidata_buffer.h>


namespace s = settings;
using namespace std;
using namespace Textra;


void class_algorithms::single_DMRG_step(class_superblock &superblock, long chi_max){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */

                    superblock.update_bond_dimensions();
    t_eig.tic();    superblock.find_ground_state(s::precision::eigSteps, s::precision::eigThreshold);    t_eig.toc();
    t_svd.tic();    superblock.truncate         (chi_max,                s::precision::SVDThreshold);    t_svd.toc();
    t_mps.tic();    superblock.update_MPS();                                                             t_mps.toc();

}

void class_algorithms::single_TEBD_step(class_superblock &superblock, long chi_max){
/*!
 * \fn single_iTEBD_step(class_superblock &superblock)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 */
                    superblock.update_bond_dimensions();
    t_evo.tic();    superblock.time_evolve();                                       t_evo.toc();
    t_svd.tic();    superblock.truncate   (chi_max, s::precision::SVDThreshold);    t_svd.toc();
    t_mps.tic();    superblock.update_MPS();                                        t_mps.toc();
}

void class_algorithms::iDMRG(){
/*!
 * \fn void iDMRG(class_superblock &superblock, class_storage &S, int max_length)
 * \brief Infinite DMRG, grows the chain from 2 up to `max_idmrg::length` particles.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param max_length Maximum chain length after which the algorithm stops.
 */

    if(!settings::idmrg::on){return;}
    ccout(0) << "Starting normal simulation using infinite-DMRG " << std::endl;
    using namespace settings::profiling;
    t_tot.set_properties(on, precision, "iDMRG Total Time           ");
    t_eig.set_properties(on, precision, "iDMRG Eigenvalue solver    ");
    t_svd.set_properties(on, precision, "iDMRG SVD Truncation       ");
    t_env.set_properties(on, precision, "iDMRG Enlarge environment  ");
    t_sto.set_properties(on, precision, "iDMRG Store MPS            ");
    t_mps.set_properties(on, precision, "iDMRG Update MPS           ");

    t_tot.tic();
    class_superblock superblock;
    class_observables observables (superblock, SimulationType::iDMRG);

    while(superblock.chain_length < s::idmrg::max_length){
        class_multidata_buffer container(&hdf5, "iDMRG/L", superblock.chain_length);

        single_DMRG_step(superblock, s::idmrg::chi_max);
                        container.push_back(observables);
        t_env.tic();    superblock.enlarge_environment();                       t_env.toc();
                        superblock.swap_AB();

    }
    t_tot.toc();
    observables.print_status_full();
    t_eig.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_env.print_time_w_percent();
    t_sto.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
    cout << endl;
}

void class_algorithms::fDMRG(){
/*!
 * \fn void fDMRG()
 * \brief Finite DMRG sweeps across the chain built during iDMRG.
 */
    if(!settings::fdmrg::on){return;}
    ccout(0) << "Starting normal simulation using finite-DMRG " << std::endl;

    using namespace settings::profiling;
    t_tot.set_properties(on, precision, "fDMRG Total Time           ");
    t_eig.set_properties(on, precision, "fDMRG Eigenvalue solver    ");
    t_svd.set_properties(on, precision, "fDMRG SVD Truncation       ");
    t_env.set_properties(on, precision, "fDMRG Enlarge environment  ");
    t_sto.set_properties(on, precision, "fDMRG Store MPS            ");
    t_mps.set_properties(on, precision, "fDMRG Update MPS           ");

    t_tot.tic();
    class_superblock     superblock;
    class_observables    observables (superblock, SimulationType::fDMRG);
    class_sweep_storage  storage(s::fdmrg::max_length);

    while(superblock.chain_length -2 < s::fdmrg::max_length){
                        single_DMRG_step(superblock, s::fdmrg::chi_max);
        t_sto.tic();    storage.insert(superblock);                  t_sto.toc();
        t_env.tic();    superblock.enlarge_environment();            t_env.toc();
                        superblock.swap_AB();
    }

    int direction  = 1;
    int sweep = 0;

    while(sweep < s::fdmrg::max_sweeps) {
                        storage.load(superblock);
                        single_DMRG_step(superblock, s::fdmrg::chi_max);
        t_sto.tic();    storage.overwrite_MPS(superblock);                                  t_sto.toc(); //Needs to occurr after update_MPS...
        t_env.tic();    superblock.enlarge_environment(direction);                          t_env.toc();
        t_sto.tic();    storage.move(superblock, direction, sweep);                         t_sto.toc();
    }

    t_tot.toc();
    observables.print_status_full();
    t_eig.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_env.print_time_w_percent();
    t_sto.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
    cout << endl;
}

void class_algorithms::iTEBD(){
/*!
 * \fn iTEBD(class_superblock &superblock, class_hdf5 &hdf5)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param max_iter Maximum number of iterations.
 */
    if(!settings::itebd::on){return;}
    ccout(0) << "Starting normal simulation using iTEBD " << std::endl;
    using namespace settings::profiling;
    class_tic_toc t_evo(on, precision, "iTEBD Time evolution       ");
    class_tic_toc t_svd(on, precision, "iTEBD SVD Truncation       ");
    class_tic_toc t_mps(on, precision, "iTEBD Update MPS           ");
    class_tic_toc t_tot(on, precision, "iTEBD Total Time           ");

    class_superblock superblock;
    class_observables observables (superblock, SimulationType::iTEBD);
    superblock.H.update_timestep(settings::itebd::delta_t, 1);
    t_tot.tic();
    for(auto step = 0; step < s::itebd::max_steps ; step++){
        class_multidata_buffer container(&hdf5, "iTEBD/step", step);
        single_TEBD_step(superblock, s::itebd::chi_max);
        if (Math::mod(step,500) == 0) { observables.print_status_update(step);}
        container.push_back(observables);
        superblock.swap_AB();

    }
    t_tot.toc();
    observables.print_status_full();
    t_evo.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
    cout << endl;
}

void class_algorithms::FES_iTEBD(){
/*!
 * \fn FES_iTEBD()
 * \brief Finite-entanglement scaling.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 *
 * This function uses infinite algorithms (iTEBD or iDMRG) until the entropy converges,
 * to see how entanglement grows as a function of  \f$\chi\f$ (chi).
 */
    if(!settings::fes_itebd::on){return;}
    using namespace settings::profiling;
    ccout(0) << "Starting Finite-Entanglement Scaling simulation using iTEBD " << std::endl;
    auto chi_max_list = Math::LinSpaced(s::fes_itebd::chi_num, s::fes_itebd::chi_min,s::fes_itebd::chi_max);
    t_tot.set_properties(on, precision, "FES_iTEBD +Total Time              ");
    t_sto.set_properties(on, precision, "          \u21B3 Store Results          ");
    t_upd.set_properties(on, precision, "          \u21B3 Update Evolution       ");
    t_sim.set_properties(on, precision, "          \u21B3+Simulation             ");
    t_evo.set_properties(on, precision, "           \u21B3 Time Evolution        ");
    t_svd.set_properties(on, precision, "           \u21B3 SVD Truncation        ");
    t_mps.set_properties(on, precision, "           \u21B3 Update MPS            ");
    for(auto &chi_max : chi_max_list) {
        t_tot.reset();
        t_tot.tic();
        t_sto.reset();
        t_upd.reset();
        t_sim.reset();
        t_evo.reset();
        t_svd.reset();
        t_mps.reset();

        double time_step = s::fes_itebd::delta_t;
        class_superblock superblock;
        class_observables observables (superblock, SimulationType::FES_iTEBD);
        class_multidata_buffer container(&hdf5, "FES_iTEBD/chi", chi_max );
        superblock.H.update_timestep(time_step,1);

        double phys_time = 0;
        int step = 0;
        double old_entropy = 0;
        double new_entropy = observables.get_entropy();
        while(time_step > 1e-6 && step < s::fes_itebd::max_steps){
            t_sim.tic();
            single_TEBD_step(superblock, chi_max );
            t_sim.toc();


            t_sto.tic();
            if(Math::mod(step,100)==0){
                container.push_back(observables);
                container.push_back(step, time_step, phys_time += time_step, t_tot.get_age());
            }
            t_sto.toc();
            if(Math::mod(step,500)==0){
//                ccout(1) << observables.first_moment() << endl;
                observables.print_status_update(step);
            }


            t_upd.tic();
            if(Math::mod(step,500)==0){
                old_entropy = new_entropy;
                new_entropy = observables.get_entropy();
                if (std::fabs((new_entropy-old_entropy)/new_entropy) < 1e-6*time_step ){
                    container.push_back(observables);
                    container.push_back(step, time_step, phys_time += time_step, t_tot.get_age());
                    time_step *=0.5;
                    superblock.H.update_timestep(time_step ,1);
                }
            }
            t_upd.toc();
            step++;
            superblock.swap_AB();

        }
        cout <<setprecision(16);
        container.push_back(observables);
        container.push_back(step, time_step, phys_time += time_step, t_tot.get_age());
        observables.print_status_full();
        t_tot.toc();
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_upd.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        t_evo.print_time_w_percent(t_sim);
        t_svd.print_time_w_percent(t_sim);
        t_mps.print_time_w_percent(t_sim);
        cout << endl;
    }
    }

void class_algorithms::FES_iDMRG(){
/*!
 * \fn FES_iDRMG()
 * \brief Finite-entanglement scaling.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 *
 * This function uses infinite algorithms (iTEBD or iDMRG) until the entropy converges,
 * to see how entanglement grows as a function of  \f$\chi\f$ (chi).
 */
    if(!settings::fes_idmrg::on){return;}
    ccout(0) << "Starting Finite-Entanglement Scaling simulation using infinite-DMRG " << std::endl;
    using namespace settings::profiling;
    t_tot.set_properties(on, precision, "FES_iDMRG +Total Time              ");
    t_sto.set_properties(on, precision, "          \u21B3 Store Results          ");
    t_env.set_properties(on, precision, "          \u21B3 Update Environments    ");
    t_sim.set_properties(on, precision, "          \u21B3+Simulation             ");
    t_eig.set_properties(on, precision, "           \u21B3 Eig. decomp.          ");
    t_svd.set_properties(on, precision, "           \u21B3 SVD Truncation        ");
    t_mps.set_properties(on, precision, "           \u21B3 Update MPS            ");

    auto chi_max_list = Math::LinSpaced(s::fes_idmrg::chi_num, s::fes_idmrg::chi_min,s::fes_idmrg::chi_max);
    for(auto &chi_max : chi_max_list ) {
        t_tot.reset();
        t_tot.tic();

        t_sim.reset();
        t_eig.reset();
        t_svd.reset();
        t_mps.reset();
        t_sto.reset();

        class_multidata_buffer container(&hdf5, "FES_iDMRG/chi", chi_max );
        class_superblock superblock;
        class_observables observables (superblock, SimulationType::FES_iDMRG);
        superblock.initialize();
        int step = 0;
        while (step  < s::fes_idmrg::max_steps) {
            t_sim.tic();
            single_DMRG_step(superblock, chi_max);
            t_sim.toc();

            if(Math::mod(step,100)==0) {
                observables.get_variance();
                observables.print_status_update(step);
            }
            t_sto.tic();
            if(Math::mod(step,100)==0) {
                container.push_back(observables);
                container.push_back(step, 0, 0, t_tot.get_age());
            }
            t_sto.toc();

            t_env.tic();
            superblock.enlarge_environment();
            t_env.toc();
            superblock.swap_AB();
            step++;
        }
        container.push_back(observables);
        container.push_back(step, 0, 0, t_tot.get_age());
        observables.print_status_full();
        t_tot.toc();
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_env.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        t_eig.print_time_w_percent(t_sim);
        t_svd.print_time_w_percent(t_sim);
        t_mps.print_time_w_percent(t_sim);

        cout << endl;
    }
}




