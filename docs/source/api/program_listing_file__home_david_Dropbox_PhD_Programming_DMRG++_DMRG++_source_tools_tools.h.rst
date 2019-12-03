
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_tools.h:

Program Listing for File tools.h
================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_tools.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/tools.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-10-13.
   //
   
   #pragma once
   
   #include <io/nmspc_logger.h>
   #include <tools/finite/measure.h>
   
   namespace tools{
       inline std::shared_ptr<spdlog::logger> log;
   }
   
   
