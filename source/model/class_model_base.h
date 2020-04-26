//
// Created by david on 2018-07-04.
//

#pragma once

#include "class_model_parameters.h"
#include <any>
#include <io/nmspc_logger.h>
#include <map>
#include <memory>
#include <simulation/enums.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace h5pp {
    class File;
}

class class_model_base {
    public:
    using Scalar   = std::complex<double>;
    using TableMap = std::map<std::string, std::any>;

    protected:
    // Common parameters
    double alpha      = 0;     /*!< Damping factor [0,1] on couplings, std::pow(J_rnd + J_ptb,1-alpha)  */
    double beta       = 0;     /*!< Damping factor [0,1] on fields, std::pow(h_rnd + h_ptb,1-alpha)  */
    double e_reduced  = 0;     /*!< "Reduced" energy offset for this mpo (to make energy-reduced MPO views) */
    double psfactor   = 0;     /*!< Parity sector separation factor */
    bool   parity_sep = false; /*!< Parity sector separation on/off */

    std::shared_ptr<spdlog::logger> log;
    Eigen::array<size_t, 4>         extent4;  /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<size_t, 2>         extent2;  /*!< Extent of pauli matrices in a rank-2 tensor */
    std::optional<size_t>           position; /*!< Position on a finite chain */
    Eigen::Tensor<Scalar, 4>        mpo_internal;
    //    [[nodiscard]] std::any &      find_val(Parameters &parameters, std::string_view key) const;
    //    [[nodiscard]] const std::any &find_val(const Parameters &parameters, std::string_view key) const;
    //    template<typename T>
    //    [[nodiscard]] T &get_val(Parameters &parameters, std::string_view key) {
    //        return std::any_cast<T &>(find_val(parameters, key));
    //    }
    //    template<typename T>
    //    [[nodiscard]] const T &get_val(const Parameters &vec, std::string_view key) {
    //        return std::any_cast<const T &>(find_val(vec, key));
    //    }

    public:
    bool all_mpo_parameters_have_been_set = false;

    explicit class_model_base(size_t position_);
    const Eigen::Tensor<Scalar, 4> &MPO() const;
    bool                            isReal() const;
    bool                            hasNaN() const;
    void                            assertValidity() const;
    void                            set_position(size_t new_pos);
    size_t                          get_position() const;
    std::vector<std::string>        get_parameter_names() const;
    std::vector<std::any>           get_parameter_values() const;
    void                            set_reduced_energy(double site_energy);
    bool                            is_damped() const;
    bool                            is_reduced() const;
    double                          get_reduced_energy() const;
    void                            print_parameter_names() const;
    void                            print_parameter_values() const;

    virtual std::unique_ptr<class_model_base> clone() const                                                                = 0;
    virtual Eigen::Tensor<Scalar, 4>          MPO_reduced_view() const                                                     = 0;
    virtual Eigen::Tensor<Scalar, 4>          MPO_reduced_view(double single_site_energy) const                            = 0;
    virtual Eigen::Tensor<Scalar, 1>          get_MPO_edge_left() const                                                    = 0;
    virtual Eigen::Tensor<Scalar, 1>          get_MPO_edge_right() const                                                   = 0;
    virtual size_t                            get_spin_dimension() const                                                   = 0;
    virtual TableMap                          get_parameters() const                                                       = 0;
    virtual void                              set_parameters(TableMap &parameters)                                          = 0;
    virtual void                              set_perturbation(double coupling_ptb, double field_ptb, PerturbMode ptbMode) = 0;
    virtual void                              set_coupling_damping(double alpha)                                           = 0;
    virtual void                              set_field_damping(double beta)                                               = 0;
    virtual void                              build_mpo()                                                                  = 0;
    virtual void                              randomize_hamiltonian()                                                      = 0;
    virtual bool                              is_perturbed() const                                                         = 0;
    virtual void                              set_averages(std::vector<TableMap> all_parameters, bool reverse = false)     = 0;
    virtual Eigen::MatrixXcd single_site_hamiltonian(size_t position, size_t sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY,
                                                     std::vector<Eigen::MatrixXcd> &SZ) const                              = 0;

    virtual ~class_model_base()                                                        = default;
    virtual void write_parameters(h5pp::File &file, std::string_view table_path) const = 0;
    virtual void read_parameters(const h5pp::File &file, std::string_view table_path)  = 0;
};
