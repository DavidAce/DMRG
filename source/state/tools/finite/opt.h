//
// Created by david on 2019-03-18.
//

#ifndef tools_FINITE_OPT_H
#define tools_FINITE_OPT_H

#include <state/tools/nmspc_tools.h>
#include <general/class_tic_toc.h>
#include <iomanip>
#include <LBFGS.h>

class class_tic_toc;


namespace tools::finite::opt{
    namespace internals{
        extern Eigen::Tensor<std::complex<double>,3> subspace_optimization       (const class_finite_state & state, const class_simulation_status & sim_status, OptType optType, OptMode optMode);
        extern Eigen::Tensor<std::complex<double>,3> direct_optimization         (const class_finite_state & state, const class_simulation_status & sim_status, OptType optType);
        extern Eigen::Tensor<std::complex<double>,3> cppoptlib_optimization      (const class_finite_state & state, const class_simulation_status & sim_status);
        extern Eigen::Tensor<std::complex<double>,4> ground_state_optimization   (const class_finite_state & state, std::string ritzstring = "SR");


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
                case OptSpace::SUBSPACE    : str << "SUBSPACE";      break;
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
        inline std::unique_ptr<class_tic_toc> t_opt  =  std::make_unique<class_tic_toc>(true,5,"t_opt ");
        inline std::unique_ptr<class_tic_toc> t_eig  =  std::make_unique<class_tic_toc>(true,5,"t_eig ");
        inline std::unique_ptr<class_tic_toc> t_ham  =  std::make_unique<class_tic_toc>(true,5,"t_ham ");
        inline std::unique_ptr<class_tic_toc> t_tot  =  std::make_unique<class_tic_toc>(true,5,"t_tot ");
        inline std::unique_ptr<class_tic_toc> t_vH2v =  std::make_unique<class_tic_toc>(true,5,"t_vH2v");
        inline std::unique_ptr<class_tic_toc> t_vHv  =  std::make_unique<class_tic_toc>(true,5,"t_vHv ");
        inline std::unique_ptr<class_tic_toc> t_vH2  =  std::make_unique<class_tic_toc>(true,5,"t_vH2 ");
        inline std::unique_ptr<class_tic_toc> t_vH   =  std::make_unique<class_tic_toc>(true,5,"t_vH  ");
        inline std::unique_ptr<class_tic_toc> t_op   =  std::make_unique<class_tic_toc>(true,5,"t_op  ");


        inline LBFGSpp::LBFGSParam<double> get_params(){
            using namespace LBFGSpp;
            LBFGSpp::LBFGSParam<double> params;
            // READ HERE http://pages.mtu.edu/~msgocken/ma5630spring2003/lectures/lines/lines/node3.html
            // I think c1 corresponds to ftol, and c2 corresponds to wolfe
            params.max_iterations = 1000;
            params.max_linesearch = 80; // Default is 20.
            params.m              = 8;     // Default is 6
            params.past           = 1;     //
            params.epsilon        = 1e-2;  // Default is 1e-5.
            params.delta          = 1e-6; // Default is 0.
            params.ftol           = 1e-4;  // Default is 1e-4.
            params.wolfe          = 0.90;   // Default is 0.9
            params.min_step       = 1e-40;
            params.max_step       = 1e+40;
            params.linesearch     = LINE_SEARCH_ALGORITHM::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
            return params;
        }

        inline auto params = get_params();

        template<typename Scalar>
        struct MultiComponents{
            Eigen::Tensor<Scalar,3> envL, envR;
            Eigen::Tensor<Scalar,4> env2L, env2R;
            Eigen::Tensor<Scalar,4> mpo;
            Eigen::Tensor<Scalar,6> mpo2;
            Eigen::DSizes<long,3>   dsizes;
            explicit MultiComponents(const class_finite_state & state);
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
            base_functor(const class_finite_state & state, const class_simulation_status &sim_status);
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
                    const class_finite_state & state,
                    const class_simulation_status &sim_status,
                    const Eigen::MatrixXcd &eigvecs_,
                    const Eigen::VectorXd  &eigvals_);

            double operator()(const Eigen::VectorXd &v_double_double, Eigen::VectorXd &grad_double_double) override;
        };


        template<typename Scalar>
        class direct_functor: public base_functor{
        private:
            using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
            using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
            MultiComponents<Scalar> multiComponents;
        public:
            explicit direct_functor(const class_finite_state & state, const class_simulation_status & sim_status);
            double operator()(const Eigen::VectorXd &v_double_double, Eigen::VectorXd &grad_double_double) override;
        };





        template<typename Scalar> using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic, 1>;

        template<typename Derived>
        VectorType<typename Derived::Scalar> get_vH2 (const Eigen::MatrixBase<Derived> &v, const MultiComponents<typename Derived::Scalar> &multiComponents){
            using Scalar = typename Derived::Scalar;
            t_vH2->tic();
            size_t log2chiL  = std::log2(multiComponents.dsizes[1]);
            size_t log2chiR  = std::log2(multiComponents.dsizes[2]);
            size_t log2spin  = std::log2(multiComponents.dsizes[0]);
            Eigen::Tensor<Scalar,3> vH2;
            if (log2spin > log2chiL + log2chiR){
                if (log2chiL > log2chiR){
                    Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), multiComponents.dsizes).shuffle(Textra::array3{1,0,2});
                    vH2 =
                            theta
                            .contract(multiComponents.env2L, Textra::idx({0}, {0}))
                            .contract(multiComponents.mpo  , Textra::idx({0,3}, {2,0}))
                            .contract(multiComponents.env2R, Textra::idx({0,3}, {0,2}))
                            .contract(multiComponents.mpo  , Textra::idx({2,1,4}, {2,0,1}))
                            .shuffle(Textra::array3{2,0,1});
                }

                else{
                    Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), multiComponents.dsizes).shuffle(Textra::array3{2,0,1});
                    vH2 =
                            theta
                            .contract(multiComponents.env2R, Textra::idx({0}, {0}))
                            .contract(multiComponents.mpo  , Textra::idx({0,3}, {2,1}))
                            .contract(multiComponents.env2L, Textra::idx({0,3}, {0,2}))
                            .contract(multiComponents.mpo  , Textra::idx({2,4,1}, {2,0,1}))
                            .shuffle(Textra::array3{2,1,0});
                }

            }else{
                Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), multiComponents.dsizes).shuffle(Textra::array3{1,0,2});
                vH2 =
                    theta
                    .contract(multiComponents.env2L, Textra::idx({0}, {0}))
                    .contract(multiComponents.mpo  , Textra::idx({0,3}, {2,0}))
                    .contract(multiComponents.mpo  , Textra::idx({4,2}, {2,0}))
                    .contract(multiComponents.env2R, Textra::idx({0,2,3}, {0,2,3}))
                    .shuffle(Textra::array3{1, 0, 2});
            }
            t_vH2->toc();
            return Eigen::Map<VectorType<Scalar>>(vH2.data(),vH2.size());
        }
        template<typename Derived>
        VectorType<typename Derived::Scalar> get_vH (const Eigen::MatrixBase<Derived> &v, const MultiComponents<typename Derived::Scalar> &multiComponents){
            using Scalar = typename Derived::Scalar;
            Eigen::Tensor<Scalar,3> vH;
            t_vH->tic();
            size_t log2chiL  = std::log2(multiComponents.dsizes[1]);
            size_t log2chiR  = std::log2(multiComponents.dsizes[2]);
//            size_t log2spin  = std::log2(multiComponents.dsizes[0]);
            if (log2chiL > log2chiR){
                Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), multiComponents.dsizes).shuffle(Textra::array3{1,0,2});
                vH =
                        theta
                        .contract(multiComponents.envL, Textra::idx({0}, {0}))
                        .contract(multiComponents.mpo , Textra::idx({0,3}, {2,0}))
                        .contract(multiComponents.envR, Textra::idx({0,2}, {0, 2}))
                        .shuffle(Textra::array3{1, 0, 2});
            }else{
                Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), multiComponents.dsizes).shuffle(Textra::array3{2,0,1});
                vH =
                        theta
                        .contract(multiComponents.envR, Textra::idx({0}, {0}))
                        .contract(multiComponents.mpo , Textra::idx({0,3}, {2,1}))
                        .contract(multiComponents.envL, Textra::idx({0,2}, {0,2}))
                        .shuffle(Textra::array3{1, 2, 0});
            }

            t_vH->toc();

            return Eigen::Map<VectorType<Scalar>>(vH.data(),vH.size());
        }

        template<typename Derived>
        std::pair<VectorType<typename Derived::Scalar>,typename Derived::Scalar>
                get_Hv_vHv(const Eigen::MatrixBase<Derived> &v,
                           const MultiComponents<typename Derived::Scalar> &multiComponents){
            auto Hv = tools::finite::opt::internals::get_vH(v,multiComponents);
            t_vHv->tic();
            auto vHv = v.dot(Hv);
            t_vHv->toc();
            return std::make_pair(Hv,vHv);
        }


        template<typename Derived>
        std::pair<VectorType<typename Derived::Scalar>,typename Derived::Scalar>
                get_H2v_vH2v(const Eigen::MatrixBase<Derived> &v,
                             const MultiComponents<typename Derived::Scalar> &multiComponents){
            auto H2v = tools::finite::opt::internals::get_vH2(v,multiComponents);
            t_vH2v->tic();
            auto vH2v = v.dot(H2v);
            t_vH2v->toc();
            return std::make_pair(H2v,vH2v);
        }
    }
}



#endif //DMRG_OPT_H
