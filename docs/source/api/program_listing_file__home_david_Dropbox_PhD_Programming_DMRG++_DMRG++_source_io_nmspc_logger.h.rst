
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_nmspc_logger.h:

Program Listing for File nmspc_logger.h
=======================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_nmspc_logger.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/io/nmspc_logger.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-03-27.
   //
   
   #pragma once
   #include <spdlog/logger.h>
   #if defined(SPDLOG_FMT_EXTERNAL)
   #include <fmt/ranges.h>
   #else
   #include <spdlog/fmt/bundled/ranges.h>
   #endif
   
   namespace Logger{
       extern void enableTimeStamp(std::shared_ptr<spdlog::logger> &log);
       extern void disableTimeStamp(std::shared_ptr<spdlog::logger> &log);
       extern void setLogLevel(std::shared_ptr<spdlog::logger> &log, size_t levelZeroToSix);
       extern size_t getLogLevel(std::shared_ptr<spdlog::logger> &log);
       extern void setLogger(std::shared_ptr<spdlog::logger> &log, std::string name, size_t levelZeroToSix = 2, bool timestamp = true);
       extern std::shared_ptr<spdlog::logger>  setLogger(std::string name, size_t levelZeroToSix = 2, bool timestamp = true);
   
   }
   
