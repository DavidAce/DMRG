
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_enums.h:

Program Listing for File enums.h
================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_enums.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/simulation/enums.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-10-16.
   //
   
   #pragma once
   enum class SimulationType      {iDMRG,fDMRG, xDMRG, iTEBD};
   enum class StorageLevel:size_t {NONE,LIGHT,NORMAL,FULL};
