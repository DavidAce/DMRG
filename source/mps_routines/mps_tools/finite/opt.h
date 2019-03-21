//
// Created by david on 2019-03-18.
//

#ifndef MPS_TOOLS_FINITE_OPT_H
#define MPS_TOOLS_FINITE_OPT_H

#include <mps_routines/nmspc_mps_tools.h>
class class_tic_toc;
class class_superblock;


namespace MPS_Tools::Finite::Opt{



    namespace internals{

        inline std::ostream& operator<<(std::ostream& str, OptMode const& mode) {
            switch (mode){
                case OptMode::OVERLAP   : str << "OVERLAP";  break;
                case OptMode::VARIANCE  : str << "VARIANCE"; break;
            }
            return str;
        }

        inline std::ostream& operator<<(std::ostream& str, OptSpace const& space) {
            switch (space){
                case OptSpace::PARTIAL  : str << "PARTIAL"; break;
                case OptSpace::FULL     : str << "FULL";    break;
                case OptSpace::DIRECT   : str << "DIRECT";  break;
            }
            return str;
        }

        extern void initialize_timers();
        extern std::shared_ptr<class_tic_toc> t_opt;
        extern std::shared_ptr<class_tic_toc> t_eig;
        extern std::shared_ptr<class_tic_toc> t_ham;
        extern std::shared_ptr<class_tic_toc> t_tot;


        struct superblock_components{
            Eigen::Tensor<double,4> HA_MPO;
            Eigen::Tensor<double,4> HB_MPO;
            Eigen::Tensor<double,3> Lblock;
            Eigen::Tensor<double,3> Rblock;
            Eigen::Tensor<double,4> Lblock2;
            Eigen::Tensor<double,4> Rblock2;
            Eigen::Tensor<double,6> HAHB;
            Eigen::Tensor<double,8> HAHB2;
            Eigen::DSizes<long,4>   dsizes;
        };



        extern double          get_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock);
        extern double          get_vHv (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock);
        extern Eigen::VectorXd get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock);
        extern Eigen::VectorXd get_vH  (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock);
        template <typename T>
        int sgn(const T val) {
            return (T(0) < val) - (val < T(0));
        }

        class base_functor{
        protected:
            double variance;
            double energy  ;
            double energy_lower_bound;
            double energy_upper_bound;
            double energy_target;
            double energy_window;
            int    counter = 0;
            bool   have_bounds_on_energy = false;

        public:
            base_functor() = default;
            double get_variance() const ;
            double get_energy  () const ;
            size_t get_count   () const ;
            void set_energy_bounds(double E_lower, double E_upper);
            virtual double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) = 0;
        };


        class subspace_functor : public base_functor {
        private:
            Eigen::MatrixXd H2;
            const Eigen::MatrixXd &eigvecs;
            const Eigen::VectorXd &eigvals;

        public:

            subspace_functor(
                    const class_superblock & superblock,
                    const Eigen::MatrixXd &eigvecs_,
                    const Eigen::VectorXd &eigvals_);

            double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) override;
        };



        class direct_functor: public base_functor{
        private:
            superblock_components superblock;
        public:
            direct_functor(const class_superblock &superblock);
            double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) override;

        };



    }





}



#endif //DMRG_OPT_H
