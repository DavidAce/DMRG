
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_mpo.cpp:

Program Listing for File mpo.cpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_mpo.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/mpo.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-06-25.
   //
   #include <state/class_state_finite.h>
   #include <tools/nmspc_tools.h>
   #include <model/class_model_base.h>
   #include <model/class_model_factory.h>
   #include <math/nmspc_random.h>
   
   void tools::finite::mpo::initialize(class_state_finite & state, const size_t length, std::string model_type){
       tools::log->trace("Initializing mpo");
       //Generate MPO
       size_t pos = 0;
       state.MPO_L.emplace_back(class_model_factory::create_mpo(pos++,model_type));
       while(true){
           state.MPO_R.emplace_back(class_model_factory::create_mpo(pos++,model_type));
           if(state.MPO_L.size() + state.MPO_R.size() >= length){break;}
       }
   }
   
   
   void tools::finite::mpo::randomize(class_state_finite &state, int seed_model) {
       tools::log->trace("Setting random fields in MPO's");
       if (seed_model >= 0){
           rn::seed(seed_model);
       }
       std::vector<std::vector<double>> all_params;
   
   
   
       for (auto &mpo : state.MPO_L){
           mpo->randomize_hamiltonian();
           all_params.push_back(mpo->get_parameter_values());
       }
       for (auto &mpo : state.MPO_R){
           mpo->randomize_hamiltonian();
           all_params.push_back(mpo->get_parameter_values());
       }
   
       //
   
       for (auto &mpo : state.MPO_L){
           mpo->set_full_lattice_parameters(all_params,false);
       }
       for (auto &mpo : state.MPO_R){
           mpo->set_full_lattice_parameters(all_params,false);
       }
   }
   
   
   
   void tools::finite::mpo::reduce_mpo_energy(class_state_finite &state){
       state.clear_measurements();
       state.clear_cache();
       if (state.active_sites.empty())
           reduce_mpo_energy_2site(state);
       else
           reduce_mpo_energy_multi(state);
   //    reduce_mpo_energy(state, state.get_multitheta());
   //
   }
   
   void tools::finite::mpo::reduce_mpo_energy_multi(class_state_finite &state){
       auto     theta = state.get_multitheta();
       double energy_per_site_before                       = tools::finite::measure::energy_per_site(state,theta);
       double energy_per_site_reduced_before               = state.get_energy_per_site_reduced();
       double energy_per_site_minus_reduced_before         = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
       double energy_variance_per_site_before              = tools::finite::measure::energy_variance_per_site(state,theta);
       log->debug("Variance check before reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(state,theta)));
       log->debug("Status before reduce (multi)         : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
                  energy_per_site_before,
                  energy_per_site_reduced_before,
                  energy_per_site_minus_reduced_before,
                  std::log10(energy_variance_per_site_before));
   
       tools::log->trace("Reducing MPO energy by: {:<20.16f}",energy_per_site_before);
       state.set_reduced_energy_per_site(energy_per_site_before);
       double energy_per_site_after                        = tools::finite::measure::energy_per_site(state,theta);
       double energy_per_site_reduced_after                = state.get_energy_per_site_reduced();
       double energy_per_site_minus_reduced_after          = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
       double energy_variance_per_site_after               = tools::finite::measure::energy_variance_per_site(state,theta);
           log->debug("Variance check after reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(state,theta)));
       log->debug("Status after reduce (multi)           : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
                  energy_per_site_after,
                  energy_per_site_reduced_after,
                  energy_per_site_minus_reduced_after,
                  std::log10(energy_variance_per_site_after));
   
       if (std::abs(energy_per_site_before - energy_per_site_after) > 0.1  )
           throw std::logic_error("Energy before and after mpo reduction differ");
       if (std::abs(energy_per_site_minus_reduced_after) > 0.1  )
           throw std::logic_error("Energy reduction failed");
   
   
   }
   
   
   void tools::finite::mpo::reduce_mpo_energy_2site(class_state_finite &state){
       state.clear_measurements();
       state.clear_cache();
       auto   theta = state.get_theta();
       double energy_per_site_before                       = tools::finite::measure::energy_per_site(state,theta);
       double energy_per_site_reduced_before               = state.get_energy_per_site_reduced();
       double energy_per_site_minus_reduced_before         = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
       double energy_variance_per_site_before              = tools::finite::measure::energy_variance_per_site(state,theta);
   //    log->debug("Variance check before reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(state)));
       log->debug("Status before reduce (2site)          : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
                  energy_per_site_before,
                  energy_per_site_reduced_before,
                  energy_per_site_minus_reduced_before,
                  std::log10(energy_variance_per_site_before));
   
   
       tools::log->trace("Reducing MPO energy by: {}",energy_per_site_before);
       state.set_reduced_energy_per_site(energy_per_site_before);
       state.clear_measurements();
       state.clear_cache();
       theta = state.get_theta();
       double energy_per_site_after                        = tools::finite::measure::energy_per_site(state,theta);
       double energy_per_site_reduced_after                = state.get_energy_per_site_reduced();
       double energy_per_site_minus_reduced_after          = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
       double energy_variance_per_site_after               = tools::finite::measure::energy_variance_per_site(state,theta);
       //    log->debug("Variance check after reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(state)));
       log->debug("Status after reduce (2site)           : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
                  energy_per_site_after,
                  energy_per_site_reduced_after,
                  energy_per_site_minus_reduced_after,
                  std::log10(energy_variance_per_site_after));
       if (std::abs(energy_per_site_before - energy_per_site_after) > 0.1  )
           throw std::logic_error("Energy before and after differ");
       if (std::abs(energy_per_site_minus_reduced_after) > 0.1  )
           throw std::logic_error("Energy reduction failed");
   }
