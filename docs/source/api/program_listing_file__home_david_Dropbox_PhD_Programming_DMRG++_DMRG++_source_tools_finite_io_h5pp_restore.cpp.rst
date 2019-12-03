
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_io_h5pp_restore.cpp:

Program Listing for File h5pp_restore.cpp
=========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_io_h5pp_restore.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/io/h5pp_restore.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-11-07.
   //
   
   #include <tools/nmspc_tools.h>
   #include <simulation/class_simulation_status.h>
   #include <state/class_state_finite.h>
   #include <model/class_model_factory.h>
   #include <h5pp/h5pp.h>
   
   void tools::finite::io::h5restore::load_from_hdf5(const h5pp::File & h5ppFile, class_state_finite & state, class_simulation_status &sim_status, const std::string & prefix_path){
       // Load into state
       try{
           sim_status = tools::common::io::h5restore::load_sim_status_from_hdf5(h5ppFile,prefix_path);
           state      = tools::finite::io::h5restore::load_state_from_hdf5(h5ppFile,prefix_path);
           state.set_sweeps(sim_status.iteration);
           tools::finite::debug::check_integrity(state);
       }catch(std::exception &ex){
           throw std::runtime_error("Failed to load from hdf5 file: " + std::string(ex.what()));
       }
   }
   
   class_state_finite tools::finite::io::h5restore::load_state_from_hdf5(const h5pp::File & h5ppFile, const std::string & prefix_path){
       class_state_finite state;
       size_t position = 0;
       size_t sites   = 0;
       Eigen::Tensor<Scalar,3> G;
       Eigen::Tensor<Scalar,1> L;
       Eigen::Tensor<Scalar,4> H;
       Eigen::MatrixXd Hamiltonian_params;
       std::string model_type;
       try{
           h5ppFile.readDataset(position             , prefix_path + "/state/position");
           h5ppFile.readDataset(sites                , prefix_path + "/state/sites");
           h5ppFile.readDataset(Hamiltonian_params   , prefix_path + "/model/Hamiltonian");
           h5ppFile.readDataset(model_type           , prefix_path + "/model/model_type");
       }catch (std::exception &ex){
           throw std::runtime_error("Couldn't read necessary model parameters: " + std::string(ex.what()));
       }
   
       try {
           for(size_t i = 0; i < sites; i++){
               h5ppFile.readDataset(G, prefix_path + "/state/mps/G_" + std::to_string(i));
               h5ppFile.readDataset(L, prefix_path + "/state/mps/L_" + std::to_string(i));
               h5ppFile.readDataset(H, prefix_path + "/state/mpo/H_" + std::to_string(i));
               if(i <= (size_t)position ) {
                   if(not state.MPS_L.empty() and state.MPS_L.back().get_chiR() != G.dimension(1)){
                       throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                   }
                   state.MPS_L.emplace_back(G,L,i);
                   state.MPO_L.emplace_back(class_model_factory::create_mpo(i,model_type,Hamiltonian_params.row(i)));
               }
               else{
                   if(not state.MPS_R.empty() and state.MPS_R.back().get_chiR() != G.dimension(1)){
                       throw std::runtime_error("Mismatch in adjacent MPS dimensions");
                   }
                   state.MPS_R.emplace_back(G,L,i);
                   state.MPO_R.emplace_back(class_model_factory::create_mpo(i,model_type,Hamiltonian_params.row(i)));
               }
           }
           h5ppFile.readDataset(state.midchain_bond()    , prefix_path + "/state/mps/L_C");
           if (state.MPS_L.size() + state.MPS_R.size() != (size_t)sites){
               throw std::runtime_error("Number of sites loaded does not match the number of sites advertised by the output file");
           }
           if (position != state.get_position()){
               throw std::runtime_error("Position loaded does not match the position read from the output file");
           }
   
       }catch (std::exception &ex){
           throw std::runtime_error("Could not read MPS/MPO tensors from file: " + std::string(ex.what()));
       }
       tools::finite::mps::rebuild_environments(state);
       return state;
   }
