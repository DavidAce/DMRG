
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_algorithm_launcher.h:

Program Listing for File class_algorithm_launcher.h
===================================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_algorithm_launcher.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/algorithms/class_algorithm_launcher.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 7/30/17.
   //
   
   #ifndef DMRG_CLASS_ALGORITHMS_H
   #define DMRG_CLASS_ALGORITHMS_H
   
   #include <memory>
   #include <spdlog/spdlog.h>
   namespace h5pp{class File;}
   
   class class_algorithm_launcher  {
   private:
       std::shared_ptr<spdlog::logger> log;
       void setLogger(std::string name);
   public:
   
   //    std::shared_ptr <class_hdf5_file> hdf5;
       std::shared_ptr<h5pp::File> h5ppFile;
       std::string hdf5_path;
       class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_);
       class_algorithm_launcher();
   
       void run_algorithms();
       void run_iDMRG();
       void run_fDMRG();
       void run_xDMRG();
       void run_iTEBD();
   
   
   };
   
   
   #endif //DMRG_CLASS_ALGORITHMS_H
