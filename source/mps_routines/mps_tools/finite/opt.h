//
// Created by david on 2019-03-18.
//

#ifndef MPS_TOOLS_FINITE_OPT_H
#define MPS_TOOLS_FINITE_OPT_H

#include <mps_routines/nmspc_mps_tools.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/meta.h>
#include <general/class_tic_toc.h>
#include <iomanip>

class class_tic_toc;
class class_superblock;
namespace LBFGSpp{
    template <typename T> class LBFGSParam;
}

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
                case OptSpace::PARTIAL     : str << "PARTIAL";      break;
                case OptSpace::FULL        : str << "FULL";         break;
                case OptSpace::DIRECT      : str << "DIRECT";       break;
                case OptSpace::GUIDED      : str << "GUIDED";       break;
                case OptSpace::CPPOPTLIB   : str << "CPPOPTLIB";    break;
            }
            return str;
        }


        namespace reports{
            using direct_opt_tuple = std::tuple<std::string,int,double,std::complex<double>,double,int,int,double>;
            using subspc_opt_tuple = std::tuple<std::string,int,double,double,double,int,int,double>;
            using lbfgs_tuple = std::tuple<double,double,double,double,double>;
            using eig_tuple   = std::tuple<int,double,double,double,double,double,double>;
//            std::vector<log_tuple> opt_log;
            void print_report(const std::vector<direct_opt_tuple> &opt_log);
            void print_report(const std::vector<subspc_opt_tuple> &opt_log);
            void print_report(const std::vector<eig_tuple> &eig_log);
            void print_report(const lbfgs_tuple lbfgs_log);
        }


        void reset_timers();
        inline std::shared_ptr<class_tic_toc> t_opt  =  std::make_shared<class_tic_toc>(true,5,"t_opt ");;
        inline std::shared_ptr<class_tic_toc> t_eig  =  std::make_shared<class_tic_toc>(true,5,"t_eig ");;
        inline std::shared_ptr<class_tic_toc> t_ham  =  std::make_shared<class_tic_toc>(true,5,"t_ham ");;
        inline std::shared_ptr<class_tic_toc> t_tot  =  std::make_shared<class_tic_toc>(true,5,"t_tot ");;
        inline std::shared_ptr<class_tic_toc> t_vH2v =  std::make_shared<class_tic_toc>(true,5,"t_vH2v");;
        inline std::shared_ptr<class_tic_toc> t_vHv  =  std::make_shared<class_tic_toc>(true,5,"t_vHv ");;
        inline std::shared_ptr<class_tic_toc> t_vH2  =  std::make_shared<class_tic_toc>(true,5,"t_vH2 ");;
        inline std::shared_ptr<class_tic_toc> t_vH   =  std::make_shared<class_tic_toc>(true,5,"t_vH  ");;
        inline std::shared_ptr<class_tic_toc> t_op   =  std::make_shared<class_tic_toc>(true,5,"t_op  ");;


        void initialize_params();
        extern std::shared_ptr<LBFGSpp::LBFGSParam<double>> params;


        struct superblock_components{
            Eigen::Tensor<double,4> HA_MPO;
            Eigen::Tensor<double,4> HB_MPO;
            Eigen::Tensor<double,3> Lblock;
            Eigen::Tensor<double,3> Rblock;
            Eigen::Tensor<double,4> Lblock2;
            Eigen::Tensor<double,4> Rblock2;
            Eigen::Tensor<double,6> HAHB;
            Eigen::Tensor<double,6> HAHA;
            Eigen::Tensor<double,6> HBHB;
            Eigen::Tensor<double,8> HAHB2;
            Eigen::Tensor<double,6> Lblock2HAHA;
            Eigen::Tensor<double,6> Rblock2HBHB;
            Eigen::DSizes<long,4>   dsizes;
        } ;



        std::pair<Eigen::VectorXd,double>    get_vH_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents);
        std::pair<Eigen::VectorXd,double>    get_vH2_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents);
        Eigen::VectorXd get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents);
        Eigen::VectorXd get_vH  (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components &superComponents);




        class base_functor{
        protected:
            double variance;
            double energy  ;
            double energy_lower_bound;
            double energy_upper_bound;
            double energy_target;
            double energy_min;
            double energy_max;
            double energy_dens;
            double energy_target_dens;
            double energy_window;
            double energy_offset;
            double norm_offset;
            size_t length;
            int    iteration;
            int    counter = 0;
            bool   have_bounds_on_energy = false;
            template <typename T>  int sgn(const T val) {return (T(0) < val) - (val < T(0)); }

            superblock_components superComponents;




        public:
            base_functor(const class_superblock & superblock, class_simulation_state &sim_state);
            double get_variance() const ;
            double get_energy  () const ;
            size_t get_count   () const ;
            virtual double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) = 0;
        };




        class direct_functor: public base_functor{
        private:
        public:
            explicit direct_functor(const class_superblock & superblock_, class_simulation_state &sim_state);
            double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) override;
        };



        class subspace_functor : public base_functor {
        private:
            Eigen::MatrixXd H2;
            const Eigen::MatrixXd &eigvecs;
            const Eigen::VectorXd &eigvals;

        public:

            explicit subspace_functor(
                    const class_superblock & superblock,
                    class_simulation_state &sim_state,
                    const Eigen::MatrixXd &eigvecs_,
                    const Eigen::VectorXd &eigvals_);

            double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) override;
        };



        class guided_functor: public base_functor{
        private:
            double windowed_func_abs(double x,double window);
            double windowed_grad_abs(double x,double window);
            double windowed_func_pow(double x,double window);
            double windowed_grad_pow(double x,double window);

        public:
            explicit guided_functor(const class_superblock &superblock, class_simulation_state & sim_state);
            double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) override;
        };






    class cppoptlib_functor: public cppoptlib::Problem<double>{
        protected:
            double variance;
            double energy  ;
            double energy_lower_bound;
            double energy_upper_bound;
            double energy_target;
            double energy_min;
            double energy_max;
            double energy_dens;
            double energy_target_dens;
            double energy_window;
            double energy_offset;
            double norm_offset;
            size_t length;
            int    iteration;
            int    counter = 0;
            bool   have_bounds_on_energy = false;
            template <typename T>  int sgn(const T val) {return (T(0) < val) - (val < T(0)); }


            superblock_components superComponents;

            double vH2v, vHv, vv,var;
            Eigen::VectorXd vH2, vH;
            bool exp_vals_computed = false;
            void compute_exp_vals (const Eigen::VectorXd &v);

        public:
            double get_variance() const ;
            double get_energy  () const ;
            size_t get_count   () const ;
            size_t get_iter    () const ;

            explicit cppoptlib_functor(const class_superblock &superblock, class_simulation_state & sim_state);
            double value (const Eigen::VectorXd &v);
            void gradient(const Eigen::VectorXd &v, Eigen::VectorXd &grad);
            bool callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x0);

        };
    }
}



#endif //DMRG_OPT_H
