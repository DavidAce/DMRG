
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_tf_ising.h:

Program Listing for File class_tf_ising.h
=========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_tf_ising.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/model/class_tf_ising.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-07-04.
   //
   
   #ifndef CLASS_TF_ISING_H
   #define CLASS_TF_ISING_H
   
   #include <iostream>
   #include <general/nmspc_tensor_extra.h>
   #include <iomanip>
   #include "class_model_base.h"
   
   class class_tf_ising : public class_model_base {
       using Scalar = std::complex<double>;
   private:
       int    spin_dim            = 0;           
       double J_coupling          = 0;
       double g_mag_field         = 0;
       double w_rnd_strength      = 0;           
       double r_rnd_field         = 0;                            
       double e_reduced           = 0;                            
       int    num_params          = 7;  //Number of parameters for this model excluding this one.
   
   public:
   
       class_tf_ising(size_t position_, std::string logName = "ISING");
       void set_hamiltonian(const Eigen::Tensor<Scalar,4> MPO, std::vector<double> parameters)  override;
       void set_hamiltonian(const std::vector<double> parameters)                               override;
       void set_hamiltonian(const Eigen::MatrixXd all_parameters, int position)                 override;
       void set_hamiltonian(const Eigen::VectorXd parameters)                                   override;
       void build_mpo()                                                                         override;
       void randomize_hamiltonian()                                                             override;
       Eigen::Tensor<Scalar,4> MPO_reduced_view()                                         const override;
       Eigen::Tensor<Scalar,4> MPO_reduced_view(double site_energy)                       const override;
       Eigen::MatrixXcd single_site_hamiltonian(
               int position,
               int sites,
               std::vector<Eigen::MatrixXcd> &SX,
               std::vector<Eigen::MatrixXcd> &SY,
               std::vector<Eigen::MatrixXcd> &SZ)                                          const override;
       std::shared_ptr<class_model_base> clone()                                     const override;
       void   set_reduced_energy(double site_energy)                                             override;
       size_t get_spin_dimension()                                                         const override;
       void   print_parameter_names ()                                                     const override;
       void   print_parameter_values()                                                     const override;
       std::vector<std::string> get_parameter_names()                                      const override;
       std::vector<double>      get_parameter_values()                                     const override;
       void   set_full_lattice_parameters(const std::vector<std::vector<double>> chain_parameters)  override;
   };
   
   #endif //CLASS_TF_ISING_H
