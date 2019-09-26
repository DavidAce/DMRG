
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_class_svd_wrapper.cpp:

Program Listing for File class_svd_wrapper.cpp
==============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_class_svd_wrapper.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/math/class_svd_wrapper.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-05-27.
   //
   
   
   
   //#define EIGEN_USE_MKL_ALL
   //#define EIGEN_USE_BLAS
   //#define EIGEN_USE_LAPACKE_STRICT
   
   //#define EIGEN_DONT_VECTORIZE
   //#define EIGEN_DONT_ALIGN_STATICALLY
   
   //#define EIGEN_MAX_ALIGN_BYTES 0
   //#define EIGEN_ENABLE_AVX512
   //#define EIGEN_UNALIGNED_VECTORIZE 0
   //#define EIGEN_DONT_ALIGN
   //#define EIGEN_MALLOC_ALREADY_ALIGNED 0
   
   
   #include <complex.h>
   #undef I
   
   //For svd debugging
   //#include <h5pp/h5pp.h>
   //#include <simulation/nmspc_settings.h>
   
   
   #include <Eigen/SVD>
   #include <Eigen/QR>
   #include <math/class_svd_wrapper.h>
   
   double class_SVD::get_truncation_error(){
       return truncation_error;
   }
   
   void class_SVD::setThreshold(double newThreshold) {
       SVDThreshold = newThreshold;
   }
   
   template<typename Scalar>
   std::tuple<class_SVD::MatrixType<Scalar>, class_SVD::VectorType<Scalar>,class_SVD::MatrixType<Scalar> , long>
   class_SVD::do_svd(const Scalar * mat_ptr, long rows, long cols, long rank_max){
       if (use_lapacke) return do_svd_lapacke(mat_ptr, rows,cols,rank_max);
   
       MatrixType<Scalar> mat = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows,cols);
   //    auto mat = Eigen::Map<const MatrixType<Scalar>>(mat_ptr, rows,cols);
       if (rows <= 0)              throw std::runtime_error("SVD error: rows() == 0");
       if (cols <= 0)              throw std::runtime_error("SVD error: cols() == 0");
       if (not mat.allFinite())    throw std::runtime_error("SVD error: matrix has inf's or nan's");
       if (mat.isZero(0))          throw std::runtime_error("SVD error: matrix is all zeros");
   
       //Debugging, save matrix into h5 file
   //    std::string outputFilename      = "svdmatrix_" + std::to_string(settings::model::seed_model) + ".h5";
   //    size_t      logLevel  = 2;
   //    h5pp::File file(outputFilename,h5pp::AccessMode::READWRITE, h5pp::CreateMode::TRUNCATE,logLevel);
   //
   //    file.writeDataset(mat, "svdmatrix");
   
   
   //    return do_svd_lapacke(mat_ptr, rows,cols,rank_max);
       Eigen::BDCSVD<MatrixType<Scalar>> SVD;
       SVD.setThreshold(SVDThreshold);
       SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
       long rank = std::min(SVD.rank(),rank_max);
       truncation_error = SVD.singularValues().normalized().tail(SVD.nonzeroSingularValues()-rank).squaredNorm();
   
       if (SVD.rank() <= 0
       or not SVD.matrixU().leftCols(rank).allFinite()
       or not SVD.singularValues().head(rank).allFinite()
       or not SVD.matrixV().leftCols(rank).allFinite() )
       {
           std::cerr   << "SVD error \n"
                       << "  SVDThreshold     = " << SVDThreshold << '\n'
                       << "  Truncation Error = " << truncation_error << '\n'
                       << "  Rank             = " << rank << '\n'
                       << "  U all finite     : " << std::boolalpha << SVD.matrixU().leftCols(rank).allFinite() << '\n'
                       << "  S all finite     : " << std::boolalpha << SVD.singularValues().head(rank).allFinite() << '\n'
                       << "  V all finite     : " << std::boolalpha << SVD.matrixV().leftCols(rank).allFinite() << '\n'
                       << "Trying SVD with LAPACKE instead \n";
           return do_svd_lapacke(mat_ptr, rows,cols,rank_max);
   //        throw std::runtime_error("SVD error:  Erroneous results");
       }
   
   //    auto [U_lapacke, S_lapacke,V_lapacke ,rank_lapacke] = do_svd_lapacke(mat_ptr, rows,cols,rank_max);
   //    long rank_common = std::min(rank_lapacke,rank);
   //    if ( (SVD.singularValues().head(rank_common) - S_lapacke).array().cwiseAbs().sum() > 1e-12 ){
   //        std::cerr   << "SVD Eigen - Lapacke mismatch" << std::endl;
   //        std::cerr   << std::setw(48) << "Eigen" << std::setw(48) << "Lapacke"  << std::setw(48) << "Diff" << std::endl;
   //        for(int i = 0; i < rank_common; i++){
   //            std::cerr   << std::setw(48) <<  std::setprecision(16)  << SVD.singularValues()(i)
   //                        << std::setw(48) <<  std::setprecision(16)  << S_lapacke(i)
   //                        << std::setw(48) <<  std::setprecision(16)  << std::abs(SVD.singularValues()(i)  - S_lapacke(i))
   //                        << std::endl;
   //        }
   //    }
   
   
   //    std::cout << "Singular values           : " << SVD.singularValues().transpose() << std::endl;
   //    std::cout << "Singular values after norm: " << SVD.singularValues().head(rank).normalized().transpose() << std::endl;
   //    std::cout << "Rank                      : " << rank << std::endl;
   //    std::cout << "Threshold                 : " << SVDThreshold << std::endl;
   //    std::cout << "Truncation error          : " << truncation_error << std::endl;
   
       return std::make_tuple(
               SVD.matrixU().leftCols(rank),
               SVD.singularValues().head(rank),
               SVD.matrixV().leftCols(rank).adjoint(),
               rank
               );
   }
   
   template std::tuple<class_SVD::MatrixType<double>, class_SVD::VectorType<double>,class_SVD::MatrixType<double> , long>
   class_SVD::do_svd(const double *, long, long, long);
   
   
   
   
   using cplx = std::complex<double>;
   template std::tuple<class_SVD::MatrixType<cplx>, class_SVD::VectorType<cplx>,class_SVD::MatrixType<cplx> , long>
   class_SVD::do_svd(const cplx *, long, long, long);
   
   
   
   
   
   template<typename Scalar>
   Eigen::Tensor<Scalar, 2>
   class_SVD::pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor){
       if (tensor.dimension(0) <= 0)  {throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(0)");}
       if (tensor.dimension(1) <= 0)  {throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(1)");}
       Eigen::Map<const MatrixType<Scalar>> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
       return Textra::Matrix_to_Tensor2(mat.completeOrthogonalDecomposition().pseudoInverse() );
   }
   
   
   
   template Eigen::Tensor<double, 2> class_SVD::pseudo_inverse(const Eigen::Tensor<double, 2> &tensor);
   template Eigen::Tensor<cplx, 2>   class_SVD::pseudo_inverse(const Eigen::Tensor<cplx  , 2> &tensor);
