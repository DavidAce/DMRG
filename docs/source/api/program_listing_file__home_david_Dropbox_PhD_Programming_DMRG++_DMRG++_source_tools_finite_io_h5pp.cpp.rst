
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_io_h5pp.cpp:

Program Listing for File h5pp.cpp
=================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_io_h5pp.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/io/h5pp.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-03-09.
   //
   
   
   #include <tools/nmspc_tools.h>
   #include <state/class_state_finite.h>
   #include <state/class_state_infinite.h>
   #include <model/class_model_factory.h>
   #include <simulation/class_simulation_status.h>
   #include <simulation/nmspc_settings.h>
   #include <general/nmspc_tensor_extra.h>
   #include <general/nmspc_quantum_mechanics.h>
   #include <stdexcept>
   #include <h5pp/h5pp.h>
   
   using Scalar    = std::complex<double>;
   
   
   bool tools::finite::io::h5dset::internals::make_extendable_dataset(const std::string & prefix_path) {
       std::string log_string = "/log";  // A logged dataset is not supposed to be extended, just written once.
       return prefix_path.find(log_string) == std::string::npos;
   }
   
   
   
   
   void tools::finite::io::h5dset::write_all_state(const class_state_finite &state, h5pp::File & h5ppFile, const std::string & prefix_path) {
       tools::log->trace("Writing state to datasets in path: {}...", prefix_path + "/state");
       tools::common::profile::t_hdf.tic();
   
       if (settings::output::storage_level >= StorageLevel::LIGHT){
           write_bond_matrix(state,h5ppFile,prefix_path);
           h5ppFile.writeDataset(state.get_length()         ,prefix_path + "/state/sites");
           h5ppFile.writeDataset(state.get_position()       ,prefix_path + "/state/position");
       }
       if(settings::output::storage_level >= StorageLevel::NORMAL){
           write_bond_matrices(state,h5ppFile,prefix_path);
       }
       if(settings::output::storage_level >= StorageLevel::FULL){
           write_full_mps(state,h5ppFile,prefix_path);
           write_full_mpo(state,h5ppFile,prefix_path);
       }
   
       tools::common::profile::t_hdf.toc();
       tools::log->trace("Writing state to datasets in path: {}... OK", prefix_path + "/state");
   
   }
   
   
   void tools::finite::io::h5dset::write_bond_matrices(const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path)
   {
       bool extendable = internals::make_extendable_dataset(prefix_path);
       auto middle = (size_t) (state.get_length() / 2);
       for (size_t i = 0; i < state.get_length(); i++){
           h5ppFile.writeDataset(state.get_MPS(i).get_L(),prefix_path + "/state/mps/L_" + std::to_string(i),extendable);
           if (i == middle){
               h5ppFile.writeDataset(state.midchain_bond(), prefix_path + "/state/mps/L_C", extendable);
           }
       }
       h5ppFile.writeDataset(state.get_truncation_errors(), prefix_path + "/state/truncation_errors",extendable);
   }
   
   
   void tools::finite::io::h5dset::write_bond_matrix(const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path)
   {
       bool extendable = internals::make_extendable_dataset(prefix_path);
       auto middle = (size_t) (state.get_length() / 2);
       h5ppFile.writeDataset(state.midchain_bond(),prefix_path + "/state/mps/L_C",extendable);
       h5ppFile.writeDataset(state.get_truncation_error(middle), prefix_path + "/state/truncation_error",extendable);
   }
   
   
   void tools::finite::io::h5dset::write_full_mps(const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path)
   {
       bool extendable = internals::make_extendable_dataset(prefix_path);
       write_bond_matrices(state,h5ppFile,prefix_path);
       for (size_t i = 0; i < state.get_length(); i++){
           h5ppFile.writeDataset(state.get_MPS(i).get_M() ,prefix_path + "/state/mps/M_" + std::to_string(i),extendable);
       }
   }
   
   
   
   
   void tools::finite::io::h5dset::write_full_mpo(const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path) {
       // Write all the MPO's
       // Remember to write tensors in row-major state order because that's what output uses.
       bool extendable = internals::make_extendable_dataset(prefix_path);
       for (auto site = 0ul; site < state.get_length(); site++){
           h5ppFile.writeDataset(state.get_MPO(site).MPO(), prefix_path + "/state/mpo/H_" + std::to_string(site),extendable);
           //Write MPO properties as attributes
           auto values = state.get_MPO(site).get_parameter_values();
           auto names  = state.get_MPO(site).get_parameter_names();
           for (size_t i = 0; i < std::min(values.size(), names.size()); i++){
               h5ppFile.writeAttributeToLink(values[i], names[i],prefix_path + "/state/mpo/H_" + std::to_string(site));
           }
       }
   }
   
   void tools::finite::io::h5dset::write_model(const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path){
       // Write down the Hamiltonian metadata as a table
       // Remember to write tensors in row-major state order because that's what output uses.
       if(settings::output::storage_level == StorageLevel::NONE) return;
       h5ppFile.writeDataset(settings::model::model_type,prefix_path + "/model/model_type");
   
       Eigen::MatrixXd hamiltonian_props;
       for (auto site = 0ul ; site < state.get_length(); site++){
           auto props = state.get_MPO(site).get_parameter_values();
           Eigen::ArrayXd  temp_row  = Eigen::Map<Eigen::ArrayXd> (props.data(),props.size());
           hamiltonian_props.conservativeResize(hamiltonian_props.rows()+1, temp_row.size());
           hamiltonian_props.bottomRows(1) = temp_row.transpose();
       }
       bool extendable = internals::make_extendable_dataset(prefix_path);
       h5ppFile.writeDataset(hamiltonian_props,prefix_path + "/model/Hamiltonian",extendable);
       int col = 0;
       for (auto &name : state.MPO_L.front()->get_parameter_names()){
           std::string attr_value = name;
           std::string attr_name  = "FIELD_" + std::to_string(col) + "_NAME";
           h5ppFile.writeAttributeToLink(attr_value, attr_name,prefix_path + "/model/Hamiltonian");
           col++;
       }
   }
   
   void tools::finite::io::h5dset::write_array_measurements(const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path){
       state.do_all_measurements();
       tools::log->trace("Writing all measurements...");
       tools::common::profile::t_hdf.tic();
       h5ppFile.writeDataset(state.measurements.bond_dimensions.value()               , prefix_path + "/bond_dimensions");
       h5ppFile.writeDataset(state.measurements.entanglement_entropies.value()        , prefix_path + "/entanglement_entropies");
       h5ppFile.writeDataset(state.measurements.spin_components.value()               , prefix_path + "/spin_components");
   
   //    h5ppFile.writeDataset(state.measurements.length.value()                        , prefix_path + "/measurements/length");
   //    h5ppFile.writeDataset(state.measurements.norm.value()                          , prefix_path + "/measurements/norm");
   //    h5ppFile.writeDataset(state.measurements.bond_dimension_midchain.value()       , prefix_path + "/measurements/bond_dimension_midchain");
   //    h5ppFile.writeDataset(state.measurements.energy.value()                        , prefix_path + "/measurements/energy");
   //    h5ppFile.writeDataset(state.measurements.energy_per_site.value()               , prefix_path + "/measurements/energy_per_site");
   //    h5ppFile.writeDataset(state.measurements.energy_variance.value()           , prefix_path + "/measurements/energy_variance");
   //    h5ppFile.writeDataset(state.measurements.energy_variance_per_site.value()      , prefix_path + "/measurements/energy_variance_per_site");
   //    h5ppFile.writeDataset(state.measurements.entanglement_entropy_midchain.value() , prefix_path + "/measurements/entanglement_entropy_midchain");
   //    h5ppFile.writeDataset(state.measurements.spin_component_sx.value()             , prefix_path + "/measurements/spin_component_sx");
   //    h5ppFile.writeDataset(state.measurements.spin_component_sy.value()             , prefix_path + "/measurements/spin_component_sy");
   //    h5ppFile.writeDataset(state.measurements.spin_component_sz.value()             , prefix_path + "/measurements/spin_component_sz");
       tools::common::profile::t_hdf.toc();
       tools::log->trace("Writing all measurements... OK");
   
   }
   
   
   
   
   
