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
public:

    class_model_base(size_t position_, std::string logName = "MODEL");
    virtual ~class_model_base() = default;
    virtual std::shared_ptr<class_model_base> clone()                     const = 0;
    const Eigen::Tensor<Scalar,4> & MPO() const;
    virtual Eigen::Tensor<Scalar,4> MPO_reduced_view()                          const = 0;
    virtual Eigen::Tensor<Scalar,4> MPO_reduced_view(double single_site_energy) const = 0;
    bool isReal()const;
    Eigen::MatrixXcd        MPO_matrix_view();    /*!< Matrix representation of full 2-site Hamiltonian */
    virtual Eigen::MatrixXcd single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY, std::vector<Eigen::MatrixXcd> &SZ) const = 0;
    virtual void   set_hamiltonian(const Eigen::Tensor<Scalar,4> MPO_, std::vector<double> parameters)  = 0;
    virtual void   set_hamiltonian(const std::vector<double> parameters)                                = 0;
    virtual void   set_hamiltonian(const Eigen::MatrixXd all_parameters, int position)                  = 0;
    virtual void   set_hamiltonian(const Eigen::VectorXd parameters)                                    = 0;
    virtual void   build_mpo()                                                                          = 0;
    virtual void   randomize_hamiltonian()                                                              = 0;
    virtual void   print_parameter_names ()                                                       const = 0;
    virtual void   print_parameter_values()                                                       const = 0;
    virtual void   set_reduced_energy(double site_energy)                                               = 0;
    virtual size_t get_spin_dimension()                                                           const = 0;
    virtual std::vector<std::string> get_parameter_names()                                        const = 0;
    virtual std::vector<double>      get_parameter_values()                                       const = 0;
    virtual void   set_full_lattice_parameters(const std::vector<std::vector<double>> parameters)          = 0;
//    virtual void   write_to_hdf5_table()                           = 0;
//    virtual double get_energy_reduced()                      const = 0;
//    virtual double get_random_field()                        const = 0;
//    virtual double get_randomness_strength()                 const = 0;
    void set_position(size_t new_pos);

    size_t    get_position()                 const;
    bool      all_mpo_parameters_have_been_set = false;
};



#endif //CLASS_HAMILTONIAN_BASE_H
