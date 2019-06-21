//
// Created by david on 2019-03-18.
//

#ifndef MPS_TOOLS_FINITE_OPT_H
#define MPS_TOOLS_FINITE_OPT_H

#include <mps_tools/nmspc_mps_tools.h>
#include <general/class_tic_toc.h>
#include <iomanip>
#include <LBFGS.h>

class class_tic_toc;
class class_superblock;


namespace mpstools::finite::opt{
    namespace internals{
        extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> subspace_optimization       (const class_superblock & superblock, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace, OptType optType);
        extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> direct_optimization(
                const class_superblock &superblock, const class_simulation_state &sim_state, OptType optType);
        extern std::tuple<Eigen::Tensor<std::complex<double>,4>, double> cppoptlib_optimization      (const class_superblock & superblock, const class_simulation_state & sim_state);
        extern std::vector<int> generate_size_list(size_t shape);


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
            }
            return str;
        }
        inline std::ostream& operator<<(std::ostream& str, OptType const& type) {
            switch (type){
                case OptType::REAL        : str << "REAL";         break;
                case OptType::CPLX        : str << "CPLX";         break;
            }
            return str;
        }

        template <typename T>  int sgn(const T val) {return (T(0) < val) - (val < T(0)); }
        double windowed_func_abs(double x,double window);
        double windowed_grad_abs(double x,double window);
        double windowed_func_pow(double x,double window);
        double windowed_grad_pow(double x,double window);



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


        inline LBFGSpp::LBFGSParam<double> get_params(){
            using namespace LBFGSpp;
            LBFGSpp::LBFGSParam<double> params;
            // READ HERE http://pages.mtu.edu/~msgocken/ma5630spring2003/lectures/lines/lines/node3.html
            // I think c1 corresponds to ftol, and c2 corresponds to wolfe
            params.max_iterations = 1000;
            params.max_linesearch = 80; // Default is 20.
            params.m              = 8;     // Default is 6
            params.past           = 1;     //
            params.epsilon        = 1e-5;  // Default is 1e-5.
            params.delta          = 1e-8; // Default is 0. Trying this one instead of ftol.
            params.ftol           = 1e-4;  // Default is 1e-4. this really helped at threshold 1e-8.
            params.wolfe          = 0.90;   // Default is 0.9
            params.min_step       = 1e-40;
            params.max_step       = 1e+40;
            params.linesearch     = LINE_SEARCH_ALGORITHM::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
            return params;
        }

        inline auto params = get_params();

        template<typename Scalar>
        struct superblock_components{
            Eigen::Tensor<Scalar,4> HA_MPO;
            Eigen::Tensor<Scalar,4> HB_MPO;
            Eigen::Tensor<Scalar,3> Lblock;
            Eigen::Tensor<Scalar,3> Rblock;
            Eigen::Tensor<Scalar,4> Lblock2;
            Eigen::Tensor<Scalar,4> Rblock2;
            Eigen::Tensor<Scalar,6> HAHB;
            Eigen::Tensor<Scalar,6> HAHA;
            Eigen::Tensor<Scalar,6> HBHB;
            Eigen::Tensor<Scalar,8> HAHB2;
            Eigen::Tensor<Scalar,6> Lblock2HAHA;
            Eigen::Tensor<Scalar,6> Rblock2HBHB;
            Eigen::DSizes<long,4>   dsizes;
            explicit superblock_components(const class_superblock &superblock);

        } ;


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
            double norm;
            size_t length;
            int    iteration;
            int    counter = 0;
            bool   have_bounds_on_energy = false;
        public:
            base_functor(const class_superblock & superblock, const class_simulation_state &sim_state);
            double get_variance() const ;
            double get_energy  () const ;
            size_t get_count   () const ;
            double get_norm    () const ;
            virtual double operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) = 0;
        };


        template <typename Scalar>
        class subspace_functor : public base_functor {
        private:
            using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
            using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
            const Eigen::MatrixXcd &eigvecs;
            const Eigen::VectorXd  &eigvals;
            MatrixType             H2;
        public:

            explicit subspace_functor(
                    const class_superblock & superblock,
                    const class_simulation_state &sim_state,
                    const Eigen::MatrixXcd &eigvecs_,
                    const Eigen::VectorXd  &eigvals_);

            double operator()(const Eigen::VectorXd &v_double_double, Eigen::VectorXd &grad_double_double) override;
        };


        template<typename Scalar>
        class direct_functor: public base_functor{
        private:
            using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
            using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
            superblock_components<Scalar> superComponents;
        public:
            explicit direct_functor(const class_superblock &superblock, const class_simulation_state & sim_state);
            double operator()(const Eigen::VectorXd &v_double_double, Eigen::VectorXd &grad_double_double) override;
        };





        template<typename Scalar> using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic, 1>;

        template<typename Derived>
        VectorType<typename Derived::Scalar> get_vH2 (const Eigen::MatrixBase<Derived> &v, const superblock_components<typename Derived::Scalar> &superComponents){
            using Scalar = typename Derived::Scalar;
            Eigen::Tensor<Scalar,4> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(v.derived().data(), superComponents.dsizes).shuffle(Textra::array4{1,0,3,2});
            t_vH2->tic();
            Eigen::Tensor<Scalar,4> vH2 =
                    theta
                            .contract(superComponents.Lblock2, Textra::idx({0}, {0}))
                            .contract(superComponents.HAHB2,   Textra::idx({5,4,0,2}, {4, 0, 1, 3}))
                            .contract(superComponents.Rblock2, Textra::idx({0,2,4}, {0, 2, 3}))
                            .shuffle(Textra::array4{1, 0, 2, 3});

            t_vH2->toc();
            return Eigen::Map<VectorType<Scalar>>(vH2.data(),vH2.size());
        }
        template<typename Derived>
        VectorType<typename Derived::Scalar> get_vH (const Eigen::MatrixBase<Derived> &v, const superblock_components<typename Derived::Scalar> &superComponents){
            using Scalar = typename Derived::Scalar;
            auto theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>> (v.derived().data(), superComponents.dsizes);
            t_vH->tic();
            Eigen::Tensor<Scalar,4> vH =
                    superComponents.Lblock
                            .contract(theta,                               Textra::idx({0},{1}))
                            .contract(superComponents.HAHB,                Textra::idx({1,2,3},{0,1,4}))
                            .contract(superComponents.Rblock ,             Textra::idx({1,3},{0,2}))
                            .shuffle(Textra::array4{1,0,2,3});
            t_vH->toc();

            return Eigen::Map<VectorType<Scalar>>(vH.data(),vH.size());
        }

        template<typename Derived>
        std::pair<VectorType<typename Derived::Scalar>,typename Derived::Scalar>
                get_Hv_vHv(const Eigen::MatrixBase<Derived> &v,
                           const superblock_components<typename Derived::Scalar> &superComponents){
            auto Hv = mpstools::finite::opt::internals::get_vH(v,superComponents);
            t_vHv->tic();
            auto vHv = v.dot(Hv);
            t_vHv->toc();
            return std::make_pair(Hv,vHv);
        }


        template<typename Derived>
        std::pair<VectorType<typename Derived::Scalar>,typename Derived::Scalar>
                get_H2v_vH2v(const Eigen::MatrixBase<Derived> &v,
                             const superblock_components<typename Derived::Scalar> &superComponents){
            auto H2v = mpstools::finite::opt::internals::get_vH2(v,superComponents);
            t_vH2v->tic();
            auto vH2v = v.dot(H2v);
            t_vH2v->toc();
            return std::make_pair(H2v,vH2v);
        }
    }
}



#endif //DMRG_OPT_H
