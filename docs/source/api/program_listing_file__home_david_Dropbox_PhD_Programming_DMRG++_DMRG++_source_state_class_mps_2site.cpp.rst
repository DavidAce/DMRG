
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_mps_2site.cpp:

Program Listing for File class_mps_2site.cpp
============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_mps_2site.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/class_mps_2site.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2017-11-13.
   //
   
   
   #include <state/class_mps_2site.h>
   #include <state/class_vidal_site.h>
   #include <math/nmspc_math.h>
   #include <general/nmspc_random_numbers.h>
   #include <simulation/nmspc_settings.h>
   #include <iomanip>
   
   using namespace std;
   using namespace Textra;
   
   using Scalar = class_mps_2site::Scalar;
   
   // Function definitions
   
   class_mps_2site::class_mps_2site(){
       MPS_A = std::make_unique<class_vidal_site>();
       MPS_B = std::make_unique<class_vidal_site>();
   }
   
   class_mps_2site::class_mps_2site(const class_mps_2site &other)
       : spin_dimension(other.spin_dimension), swapped(other.swapped), MPS_A(copy_unique(other.MPS_A)), MPS_B(copy_unique(other.MPS_B)),  LC(other.LC)
   {
   }
   
   
   bool class_mps_2site::isReal()const {return MPS_A->isReal() and MPS_B->isReal();}
   long class_mps_2site::chiA () const {return MPS_A->get_L().dimension(0);}
   long class_mps_2site::chiB () const {return MPS_B->get_L().dimension(0);}
   long class_mps_2site::chiC () const {return LC.dimension(0);}
   long class_mps_2site::spindim() const {return spin_dimension;}
   Eigen::Tensor<Scalar,3> class_mps_2site::A() const {return MPS_A->get_A();}
   Eigen::Tensor<Scalar,3> class_mps_2site::B() const {return MPS_B->get_B();}
   Eigen::Tensor<Scalar,2> class_mps_2site::C() const {return Textra::asDiagonal(LC);}
   void class_mps_2site::set_mps(const Eigen::Tensor<Scalar,1> &LA,
                                 const Eigen::Tensor<Scalar,3> &GA,
                                 const Eigen::Tensor<Scalar,1> &LC_,
                                 const Eigen::Tensor<Scalar,3> &GB,
                                 const Eigen::Tensor<Scalar,1> &LB)
   {
       MPS_A->set_mps(GA,LA);
       MPS_B->set_mps(GB,LB);
       LC = LC_;
   }
   
   Eigen::DSizes<long,4> class_mps_2site::dimensions() const {return Eigen::DSizes<long,4>{spin_dimension,chiA(), spin_dimension,chiB()};}
   
   
   
   void class_mps_2site::initialize(int spin_dim){
       spin_dimension = spin_dim;
       Eigen::Tensor<Scalar,3> GA(spin_dimension,1,1);
       Eigen::Tensor<Scalar,3> GB(spin_dimension,1,1);
       Eigen::Tensor<Scalar,1> LA(1);
       Eigen::Tensor<Scalar,1> LC(1);
       Eigen::Tensor<Scalar,1> LB(1);
   
       // Default is a product state, spins pointing up in z.
       GA.setZero();
       GB.setZero();
       LA.setConstant(1.0);
       LB.setConstant(1.0);
       LC.setConstant(1.0);
   
       GA(0,0,0) = 1;
       GB(0,0,0) = 1;
   
       set_mps(LA,GA,LC,GB,LB);
       swapped = false;
   //    theta = get_theta();
   }
   
   
   void class_mps_2site::swap_AB() {
       swapped = !swapped;
   
       //Swap Gamma
       tmp3 = MPS_A->get_G();
       MPS_A->set_G(MPS_B->get_G());
       MPS_B->set_G(tmp3);
   
   
       tmp1 = LC;
       LC = MPS_B->get_L();
       MPS_A->set_L(tmp1);
       MPS_B->set_L(tmp1);
   }
   
   Eigen::Tensor<class_mps_2site::Scalar,4> class_mps_2site::get_theta(Scalar norm) const
   {
       return
               A().contract(C(), idx({2},{0}))
                  .contract(B(), idx({2},{1})) / norm;
   }
