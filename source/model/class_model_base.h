//
// Created by david on 2018-07-04.
//

#ifndef CLASS_HAMILTONIAN_BASE_H
#define CLASS_HAMILTONIAN_BASE_H

#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include <io/nmspc_logger.h>

class class_model_base{
    using Scalar = std::complex<double>;
protected:
    Eigen::array<long, 4> extent4;                           /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2;                           /*!< Extent of pauli matrices in a rank-2 tensor */
    std::optional<size_t> position;                          /*!< Position on a finite chain */
    std::shared_ptr<spdlog::logger> log;
    Eigen::Tensor<Scalar,4> mpo_internal;
    double e_reduced = 0;                                   /*!< "Reduced" energy offset for this mpo (to make "reduced" MPO views) */

public:

    explicit class_model_base(size_t position_, std::string logName = "MODEL");
    bool   all_mpo_parameters_have_been_set =   false;
    const Eigen::Tensor<Scalar,4> &  MPO()      const;
    bool   isReal()                             const;
    void   set_position(size_t new_pos);
    size_t get_position()                       const;

    Eigen::MatrixXcd MPO_matrix_view();                /*!< Matrix representation of full 2-site Hamiltonian */

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



#endif //CLASS_HAMILTONIAN_BASE_H
