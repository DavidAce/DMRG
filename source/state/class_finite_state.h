//
// Created by david on 2019-01-29.
//

#ifndef CLASS_FINITE_CHAIN_STORAGE_H
#define CLASS_FINITE_CHAIN_STORAGE_H
#include <memory>
#include <complex>
#include <optional>
#include <general/nmspc_tensor_extra.h>
#include <state/class_environment.h>
#include <state/class_vidal_site.h>
#include <model/class_model_base.h>




/**
 * \class class_finite_state
 *
 * \brief Stores the finite mps state and components to operate on it, such as environments and model mpo's.
 *
 * The finite state mps is always partitioned into two lists, corresponding to the two sides left and right of the current position.
 * Each site contains a "Vidal" site, i.e., a bond matrix \f$ \Lambda \f$ and an mps tensor \f$\Gamma \f$ in Vidal's notation.
 * In the left side the sites are left normalized, with the bond matrix on the left \f$ \Lambda \Gamma \ = A\f$.
 * In the right side the sites are right normalized, with the bond matrix on the right \f$ \Gamma \Lambda = B\f$.
 * At the center is the center bond matrix MPS_C, a.k.a. \f$ \Lambda^C \f$. For a state with 10 sites the layout is seen below
 *
 * \code
 * |Finite state> =  ---MPS_L(0)--- ... MPS_L(4)---MPS_L(5)---MPS_C---MPS_R(6)---MPS_R(7)---...---MPS_R(9)
 * \endcode
 *
 *  The numbers in parentheses denote the position in the chain, note that this isn't the same as the position in the corresponding containers.
 */

class class_finite_state
{
public:
    using Scalar = std::complex<double>;
private:
    using MType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    using VType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    template<auto rank> using TType = Eigen::Tensor<Scalar,rank>;
    int num_sweeps   = 0;
    int num_moves    = 0;
    int direction    = 1;
    long chi_max     = 0;
public:
    class_finite_state()=default;
    ~class_finite_state();
    class_finite_state(const class_finite_state & other);
    class_finite_state& operator= (const class_finite_state & other);



    std::list<class_vidal_site>                        MPS_L;   /*!< A list of stored \f$ \Lambda^B \Gamma^A...  \f$-tensors. */
    std::list<class_vidal_site>                        MPS_R;   /*!< A list of stored \f$ \Gamma^B \Lambda^B...  \f$-tensors. */
    TType<1>                                           MPS_C;   /*!< Current center bond matrix. */
    std::list<class_environment>                       ENV_L;
    std::list<class_environment>                       ENV_R;
    std::list<class_environment_var>                   ENV2_L;
    std::list<class_environment_var>                   ENV2_R;
    std::list<std::unique_ptr<class_model_base>>       MPO_L;     /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::list<std::unique_ptr<class_model_base>>       MPO_R;     /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */

    void do_all_measurements();



    int    get_sweeps() const ;
    int    reset_sweeps();
    void   set_sweeps(int num_sweeps_);
    void   increment_sweeps();

    int    get_moves() const ;
    int    reset_moves();
    void   set_moves(int num_moves_);
    void   increment_moves();



    long   get_chi_max()                        const;
    void   set_chi_max(long chi_max_);
    void set_positions();
    size_t get_length()                         const;
    size_t get_position()                       const;
    void flip_direction();
    int  get_direction()                        const;
    Eigen::DSizes<long,3>  dimensions_2site()   const;
    size_t size_2site()                         const;


    bool position_is_the_middle()               const;
    bool position_is_the_middle_any_direction() const;
    bool position_is_the_left_edge()            const;
    bool position_is_the_right_edge()           const;
    bool position_is_any_edge()                 const;
    bool position_is_at(size_t pos)             const;
    bool isReal()                               const;


    const class_vidal_site       & get_MPS(size_t pos)                  const;
          class_vidal_site       & get_MPS(size_t pos);
    const class_model_base       & get_MPO(size_t pos)                  const;
          class_model_base       & get_MPO(size_t pos);
    const class_environment      & get_ENVL(size_t pos)                 const;
    const class_environment      & get_ENVR(size_t pos)                 const;
    const class_environment_var  & get_ENV2L(size_t pos)                const;
    const class_environment_var  & get_ENV2R(size_t pos)                const;
    const Eigen::Tensor<Scalar,3> & get_G(size_t pos)                   const;
    const Eigen::Tensor<Scalar,1> & get_L(size_t pos)                   const;
    TType<3> & get_G(size_t pos);
    TType<1> & get_L(size_t pos);
    TType<3>   get_A()                                   const;
    TType<3>   get_B()                                   const;
    TType<3>   get_A(size_t pos)                         const;
    TType<3>   get_B(size_t pos)                         const;
    TType<4>   get_theta()                               const;
    TType<4>   get_theta(size_t pos)                     const;

    // For reduced energy MPO's
    bool   isReduced()                            const;
    double get_energy_reduced()                   const;
    void   set_reduced_energy(double site_energy);



    //For multisite
    std::list<size_t>      active_sites;
    std::list<size_t>      activate_sites(const long threshold, const size_t max_sites);
    Eigen::DSizes<long,3>  active_dimensions() const;
    size_t                 active_problem_size() const;

    const TType<3> & get_multitheta()  const;
    const TType<4> & get_multimpo  ()  const;
    std::pair<std::reference_wrapper<const class_environment>     , std::reference_wrapper<const class_environment>>      get_multienv ()     const;
    std::pair<std::reference_wrapper<const class_environment_var> , std::reference_wrapper<const class_environment_var>>  get_multienv2()     const;

//    TType<6>   get_multi_hamiltonian() const;
//    TType<6>   get_multi_hamiltonian2() const;
//    MType get_multi_hamiltonian_matrix() const;
//    MType get_multi_hamiltonian2_matrix() const;
//    MType get_multi_hamiltonian2_subspace_matrix(const MType & eigvecs ) const;
//

    std::vector<double>  truncation_error;
    struct Measurements {
        std::optional<size_t>               length                                  = {};
        std::optional<size_t>               bond_dimension_midchain                 = {};
        std::optional<size_t>               bond_dimension_current                  = {};
        std::optional<std::vector<size_t>>  bond_dimensions                         = {};
        std::optional<double>               norm                                    = {};
        std::optional<double>               energy                                  = {};
        std::optional<double>               energy_per_site                         = {};
        std::optional<double>               energy_variance_mpo                     = {};
        std::optional<double>               energy_variance_per_site                = {};
        std::optional<double>               spin_component_sx                       = {};
        std::optional<double>               spin_component_sy                       = {};
        std::optional<double>               spin_component_sz                       = {};
        std::optional<std::vector<double>>  spin_components                         = {};
        std::optional<double>               entanglement_entropy_midchain           = {};
        std::optional<double>               entanglement_entropy_current            = {};
        std::optional<std::vector<double>>  entanglement_entropies                  = {};
    };
    mutable Measurements measurements;
    void unset_measurements() const;
    void do_all_measurements()const;
    void clear_cache() const;

    void tag_active_sites_have_been_updated(bool tag)   const;
    void tag_all_sites_have_been_updated(bool tag)      const;
    bool all_sites_updated() const;
    bool active_sites_updated() const;
    mutable std::vector<bool> site_update_tags;
private:
    struct Cache{
        std::optional<TType<4>> multimpo        = {};
        std::optional<TType<3>> multitheta      = {};
//        std::optional<TType<6>> multiham        = {};
//        std::optional<TType<6>> multiham_sq     = {};
//        std::optional<MType>    multiham_mat    = {};
//        std::optional<MType>    multiham_sq_mat = {};
//        std::optional<MType>    multiham_sq_sub = {};
    };
    mutable Cache cache;
};


#endif //DMRG_CLASS_FINITE_CHAIN_STORAGE_H
