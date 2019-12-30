
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_state_infinite.h:

Program Listing for File class_state_infinite.h
===============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_state_infinite.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/class_state_infinite.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 7/22/17.
   //
   
   #pragma once
   
   
   #include <memory>
   #include <optional>
   #include <general/nmspc_tensor_extra.h>
   #include <general/class_tic_toc.h>
   #include <math/nmspc_eigutils.h>
   #include <simulation/nmspc_settings.h>
   #include <io/nmspc_logger.h>
   class class_mps_2site;
   class class_model_base;
   class class_environment;
   class class_environment_var;
   
   class class_state_infinite {
   public:
       SimulationType sim_type;
       std::string sim_name;
   private:
       std::optional<long> chi_lim;
       std::optional<long> chi_max;
       std::shared_ptr<spdlog::logger> log;
   public:
       using Scalar = std::complex<double>;
   
       class_state_infinite(SimulationType sim_type_, std::string sim_name_);
   
   //    // we can use the default move constructor
   //    class_state_infinite(class_state_infinite&&);
   //    // We have to provide a deep copying assignment operator
   //    class_state_infinite& operator=(class_state_infinite const& source);
   //    // we can use the default move assignment operator
   //    class_state_infinite& operator=(class_state_infinite&&);
   //    ~class_state_infinite();
       void clear();
   
   
       std::shared_ptr<class_mps_2site>         MPS;        
       std::shared_ptr<class_model_base>        HA;         
       std::shared_ptr<class_model_base>        HB;         
       std::shared_ptr<class_environment>       Lblock;     
       std::shared_ptr<class_environment>       Rblock;     
       std::shared_ptr<class_environment_var>   Lblock2;    
       std::shared_ptr<class_environment_var>   Rblock2;    
       size_t                   get_length()     const;
       size_t                   get_position()   const;
       long                     get_chi()        const ;
       long                     get_chi_lim()    const;
       void                     set_chi_lim(long chi_lim_);
       long                     get_chi_max()    const;
       void                     set_chi_max(long chi_max_);
       double                   get_truncation_error()const;
       Eigen::Tensor<Scalar, 4> get_theta()      const;
       Eigen::DSizes<long,4>    dimensions()     const;
       void assert_positions() const;
   //
   //    Eigen::Tensor<Scalar, 4>
   //    optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::Ritz::SR
   //    )    __attribute((hot));                            /*!< Finds the smallest algebraic eigenvalue and eigenvector (the ground state) using [Spectra](https://github.com/yixuan/spectra). */
   //
   //
   //    Eigen::Tensor<Scalar, 4>
   //    evolve_MPS(const Eigen::Tensor<Scalar, 4> &U);
   //    Eigen::Tensor<Scalar, 4>
   //    evolve_MPS(const Eigen::Tensor<Scalar, 4> &theta, const Eigen::Tensor<Scalar, 4> &U);
   //    Eigen::Tensor<Scalar,4> truncate_MPS(
   //            const Eigen::Tensor<Scalar, 4> &theta,        /*!< The 2-site MPS to truncate */
   //            long chi_,                                      /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
   //            double svd_threshold                            /*!< Minimum threshold value for keeping singular values. */
   //    )             __attribute((hot));                      /*!< Singular value decomposition, SVD, or Schmidt decomposition, of the ground state, where the truncation keeps \f$\chi\f$ (`chi`) singular values. */
   
   //    void truncate_MPS(
   //            const Eigen::Tensor<Scalar, 4> &theta,        /*!< The 2-site MPS to truncate */
   //            class_mps_2site &MPS_out,
   //            long chi_,                                      /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
   //            double svd_threshold                            /*!< Minimum threshold value for keeping singular values. */
   //    )             __attribute((hot));                      /*!< Singular value decomposition, SVD, or Schmidt decomposition, of the ground state, where the truncation keeps \f$\chi\f$ (`chi`) singular values. */
   
   
   
      void enlarge_environment(int direction = 0);          
       bool isReal() const;
       template<typename T>  Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> get_H_local_matrix()            const;
       template<typename T>  Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> get_H_local_sq_matrix ()        const;
   
   
   
       void set_superblock(
               const Eigen::Tensor<Scalar,4> &Lblock2_,
               const Eigen::Tensor<Scalar,3> &Lblock_,
               const Eigen::Tensor<Scalar,4> &MPO_A,
               const Eigen::Tensor<Scalar,1> &LA,
               const Eigen::Tensor<Scalar,3> &GA,
               const Eigen::Tensor<Scalar,1> &LC,
               const Eigen::Tensor<Scalar,3> &GB,
               const Eigen::Tensor<Scalar,1> &LB,
               const Eigen::Tensor<Scalar,4> &MPO_B,
               const Eigen::Tensor<Scalar,3> &Rblock_,
               const Eigen::Tensor<Scalar,4> &Rblock2_
               );
   
       void set_positions(int position);
   
   //    void set_current_dimensions()      ;                /*!< Update variables for dimensions */
       void swap_AB();                                     
       struct Measurements {
           std::optional<size_t> length                            = {};
           std::optional<size_t> bond_dimension                    = {};
           std::optional<double> current_entanglement_entropy      = {};
           std::optional<double> norm                              = {};
           std::optional<double> energy_mpo                        = {};
           std::optional<double> energy_per_site_mpo               = {};
           std::optional<double> energy_variance_mpo               = {};
           std::optional<double> energy_per_site_ham               = {};
           std::optional<double> energy_per_site_mom               = {};
           std::optional<double> energy_variance_per_site_mpo      = {};
           std::optional<double> energy_variance_per_site_ham      = {};
           std::optional<double> energy_variance_per_site_mom      = {};
           std::optional<double> truncation_error                  = {};
       };
       mutable Measurements measurements;
       mutable bool has_been_written  = false;
       void do_all_measurements() const;
       void unset_measurements()  const;
   
   
   
   
       //Profiling
   //    mutable class_tic_toc t_eig;
   //    mutable class_tic_toc t_ene_mpo;
   //    mutable class_tic_toc t_ene_ham;
   //    mutable class_tic_toc t_ene_mom;
   //    mutable class_tic_toc t_var_mpo;
   //    mutable class_tic_toc t_var_ham;
   //    mutable class_tic_toc t_var_mom;
   //    mutable class_tic_toc t_entropy;
   //    mutable class_tic_toc t_temp1;
   //    mutable class_tic_toc t_temp2;
   //    mutable class_tic_toc t_temp3;
   //    mutable class_tic_toc t_temp4;
   
   //    void set_profiling_labels();
   //    void print_profiling(class_tic_toc &t_parent);
   
   
   
   };
   
   
