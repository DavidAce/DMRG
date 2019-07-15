
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_xDMRG.h:

Program Listing for File class_xDMRG.h
======================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_xDMRG.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/algorithms/class_xDMRG.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-02-09.
   //
   
   #ifndef DMRG_CLASS_EXITED_DMRG_H
   #define DMRG_CLASS_EXITED_DMRG_H
   
   #include "class_algorithm_finite.h"
   #include <unsupported/Eigen/CXX11/Tensor>
   #include <Eigen/Core>
   
   class class_log_dmrg;
   
   
   
   
   class class_finite_state;
   class class_xDMRG : public class_algorithm_finite {
   private:
   
   public:
       //Inherit the constructor of class_algorithm_base
       using class_algorithm_finite::class_algorithm_finite;
       explicit class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_);
       std::unique_ptr<class_hdf5_log<class_log_dmrg>> log_dmrg;
   
       bool   projected_during_saturation  = false;
   
       void find_energy_range();
       void single_DMRG_step();
       void run_preprocessing()                final;
       void run_simulation()                   final;
       void check_convergence()                final;
       void write_logs(bool force = false)     final;
       bool   sim_on()                         final;
       long   chi_max()                        final;
       size_t num_sites()                      final;
       size_t write_freq()                     final;
       size_t print_freq()                     final;
       bool   chi_grow()                       final;
       bool   store_wave_function()            final;
   
   };
   
   
   
   #endif //DMRG_CLASS_EXITED_DMRG_H
