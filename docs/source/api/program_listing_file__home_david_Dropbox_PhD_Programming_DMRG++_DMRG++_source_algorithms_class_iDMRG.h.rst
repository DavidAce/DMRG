
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_iDMRG.h:

Program Listing for File class_iDMRG.h
======================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_iDMRG.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/algorithms/class_iDMRG.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-01-18.
   //
   
   #ifndef DMRG_CLASS_INFINITE_DMRG_H
   #define DMRG_CLASS_INFINITE_DMRG_H
   
   
   #include "class_algorithm_infinite.h"
   class class_log_dmrg;
   
   class class_iDMRG : public class_algorithm_infinite {
   public:
       //Inherit the constructor of class_algorithm_base
       using class_algorithm_infinite::class_algorithm_infinite;
       explicit class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_);
       std::unique_ptr<class_hdf5_log<class_log_dmrg>> log_dmrg;
   
       void single_DMRG_step(std::string ritz);
       void run_simulation()                                   final;
       void check_convergence()                                final;
       void write_logs(bool force = false)                     final;
       bool    sim_on ()                                       final;
       long    chi_max()                                       final;
       size_t  num_sites()                                     final;
       size_t  write_freq()                                    final;
       size_t  print_freq()                                    final;
       bool    chi_grow()                                      final;
   };
   
   
   #endif //DMRG_CLASS_INFINITE_DMRG_H
