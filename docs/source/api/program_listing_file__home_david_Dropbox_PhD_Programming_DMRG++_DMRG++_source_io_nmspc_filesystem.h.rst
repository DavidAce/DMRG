
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_nmspc_filesystem.h:

Program Listing for File nmspc_filesystem.h
===========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_nmspc_filesystem.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/io/nmspc_filesystem.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-12-09.
   //
   
   #ifndef DMRG_NMSPC_FILESYSTEM_H
   #define DMRG_NMSPC_FILESYSTEM_H
   
   // Include filesystem or experimental/filesystem
   #if __has_include(<filesystem>)
   #include <filesystem>
       namespace tools{
           namespace fs =  std::filesystem;
       }
   #elif __has_include(<experimental/filesystem>)
   #include <experimental/filesystem>
   namespace tools {
       namespace fs = std::experimental::filesystem;
   }
   #else
   #error Could not find <filesystem> or <experimental/filesystem>
   #endif
   
   
   #endif //DMRG_NMSPC_FILESYSTEM_H
