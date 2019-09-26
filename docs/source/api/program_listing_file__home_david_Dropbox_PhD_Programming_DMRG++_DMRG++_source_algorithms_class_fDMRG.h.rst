
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_fDMRG.h:

Program Listing for File class_fDMRG.h
======================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_fDMRG.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/algorithms/class_fDMRG.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-01-31.
   //
   
   #ifndef DMRG_CLASS_FINITE_DMRG_H
   #define DMRG_CLASS_FINITE_DMRG_H
   
   #include "class_algorithm_finite.h"
   class class_log_finite_dmrg_measurements;
   
   
   class class_finite_state;
   class class_fDMRG : public class_algorithm_finite {
   public:
       //Inherit the constructor of class_algorithm_base
       using class_algorithm_finite::class_algorithm_finite;
       explicit class_fDMRG(std::shared_ptr<h5pp::File> h5ppFile_);
       std::unique_ptr<class_hdf5_log<class_log_finite_dmrg_measurements>> log_dmrg;
   
       bool   projected_during_saturation  = false;
       void run_simulation()                                        final;
       void check_convergence()                                     final;
       bool   sim_on()                                              final;
       long   chi_max()                                             final;
       size_t num_sites()                                           final;
       size_t write_freq()                                          final;
       size_t print_freq()                                          final;
       bool   chi_grow()                                            final;
       bool   store_wave_function()                                 final;
   
   };
   
   
   
   #endif //DMRG_CLASS_FINITE_DMRG_H
