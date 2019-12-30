
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_base.h:

Program Listing for File class_model_base.h
===========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_base.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/model/class_model_base.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-07-04.
   //
   
   #pragma once
   
   #include <memory>
   #include <unsupported/Eigen/CXX11/Tensor>
   #include <io/nmspc_logger.h>
   
   class class_model_base{
       using Scalar = std::complex<double>;
   protected:
       Eigen::array<long, 4> extent4;                           
       Eigen::array<long, 2> extent2;                           
       std::optional<size_t> position;                          
       std::shared_ptr<spdlog::logger> log;
       Eigen::Tensor<Scalar,4> mpo_internal;
       double e_reduced = 0;                                   
   public:
   
       explicit class_model_base(size_t position_, std::string logName = "MODEL");
       bool   all_mpo_parameters_have_been_set =   false;
       const Eigen::Tensor<Scalar,4> &  MPO()      const;
       bool   isReal()                             const;
       void   set_position(size_t new_pos);
       size_t get_position()                       const;
   
       Eigen::MatrixXcd MPO_matrix_view();                
       void set_reduced_energy(double site_energy);
       bool isReduced ()                           const;
       double get_reduced_energy()                 const;
   
       virtual ~class_model_base() = default;
       virtual std::unique_ptr<class_model_base> clone()                                                       const = 0;
       virtual Eigen::Tensor<Scalar,4> MPO_reduced_view()                                                      const = 0;
       virtual Eigen::Tensor<Scalar,4> MPO_reduced_view(double single_site_energy)                             const = 0;
       virtual void   set_hamiltonian(const Eigen::Tensor<Scalar,4> & MPO_, std::vector<double> & parameters)        = 0;
       virtual void   set_hamiltonian(const std::vector<double> & parameters)                                        = 0;
       virtual void   set_hamiltonian(const Eigen::MatrixXd & all_parameters, int position)                          = 0;
       virtual void   set_hamiltonian(const Eigen::VectorXd & parameters)                                            = 0;
       virtual void   build_mpo()                                                                                    = 0;
       virtual void   randomize_hamiltonian()                                                                        = 0;
       virtual void   perturb_hamiltonian(double amplitude)                                                          = 0;
       virtual bool   is_perturbed()                                                                           const = 0;
       virtual void   print_parameter_names ()                                                                 const = 0;
       virtual void   print_parameter_values()                                                                 const = 0;
       virtual size_t get_spin_dimension()                                                                     const = 0;
       virtual std::vector<std::string> get_parameter_names()                                                  const = 0;
       virtual std::vector<double>      get_parameter_values()                                                 const = 0;
       virtual void
       set_full_lattice_parameters(const std::vector<std::vector<double>> parameters, bool reverse = false) = 0;
       virtual Eigen::MatrixXcd single_site_hamiltonian(
               int position,
               int sites,
               std::vector<Eigen::MatrixXcd> &SX,
               std::vector<Eigen::MatrixXcd> &SY,
               std::vector<Eigen::MatrixXcd> &SZ)                                                    const = 0;
   
   };
   
   
   
