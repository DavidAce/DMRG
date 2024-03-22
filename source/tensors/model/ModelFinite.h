#pragma once
#include "config/enums.h"
#include "math/float.h"
#include "math/svd/config.h"
#include <any>
#include <array>
#include <complex>
#include <memory>
#include <unordered_map>
#include <unsupported/Eigen/CXX11/Tensor>

class MpoSite;
class TensorsFinite;
class ModelLocal;
enum class MposWithEdges { OFF, ON };

class ModelFinite {
    private:
    friend TensorsFinite;
    struct Cache {
        std::optional<std::vector<size_t>>                        cached_sites              = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 4>>                     multisite_mpo             = std::nullopt;
        std::optional<std::vector<size_t>>                        multisite_mpo_ids         = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 2>>                     multisite_ham             = std::nullopt;
        std::optional<Eigen::Tensor<cplx_t, 4>>                   multisite_mpo_t           = std::nullopt;
        std::optional<Eigen::Tensor<cplx_t, 2>>                   multisite_ham_t           = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 4>>                     multisite_mpo_squared     = std::nullopt;
        std::optional<std::vector<size_t>>                        multisite_mpo_squared_ids = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 2>>                     multisite_ham_squared     = std::nullopt;
        std::unordered_map<std::string, Eigen::Tensor<cplx, 4>>   multisite_mpo_temps;   // Keeps previous results for reuse
        std::unordered_map<std::string, Eigen::Tensor<cplx_t, 4>> multisite_mpo_t_temps; // Keeps previous results for reuse
    };
    mutable Cache cache;
    //    std::vector<Eigen::Tensor<cplx, 4>> get_compressed_mpos(std::vector<Eigen::Tensor<cplx, 4>> mpos);
    void                             randomize();
    void                             clear_mpo_squared();
    bool                             has_mpo_squared() const;
    void                             set_energy_shift(double total_energy);
    void                             set_energy_shift_per_site(double energy_shift_per_site);
    bool                             set_parity_shift_mpo(int sign, std::string_view axis);
    std::pair<int, std::string_view> get_parity_shift_mpo() const;
    bool                             has_parity_shift_mpo() const;
    bool                             set_parity_shift_mpo_squared(int sign, std::string_view axis);
    std::pair<int, std::string_view> get_parity_shift_mpo_squared() const;
    bool                             has_parity_shift_mpo_squared() const;

    public:
    std::vector<std::unique_ptr<MpoSite>> MPO; /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::vector<size_t>                   active_sites;
    ModelType                             model_type = ModelType::ising_tf_rf;

    public:
                 ModelFinite();
    ~            ModelFinite();                         // Read comment on implementation
                 ModelFinite(ModelFinite &&other);      // default move ctor
    ModelFinite &operator=(ModelFinite &&other);        // default move assign
                 ModelFinite(const ModelFinite &other); // copy ctor
    ModelFinite &operator=(const ModelFinite &other);   // copy assign

    void           initialize(ModelType model_type_, size_t model_size);
    size_t         get_length() const;
    bool           is_real() const;
    bool           has_nan() const;
    void           assert_validity() const;
    const MpoSite &get_mpo(size_t pos) const;
    MpoSite       &get_mpo(size_t pos);
    void           build_mpo();
    void           build_mpo_squared();
    void           compress_mpo();
    void           compress_mpo_squared();

    std::vector<std::reference_wrapper<const MpoSite>> get_mpo(const std::vector<size_t> &sites) const;
    std::vector<std::reference_wrapper<MpoSite>>       get_mpo(const std::vector<size_t> &sites);
    std::vector<std::reference_wrapper<const MpoSite>> get_mpo_active() const;
    std::vector<std::reference_wrapper<MpoSite>>       get_mpo_active();
    std::vector<Eigen::Tensor<cplx, 4>>                get_all_mpo_tensors(MposWithEdges withEdges = MposWithEdges::OFF);
    std::vector<Eigen::Tensor<cplx_t, 4>>              get_all_mpo_tensors_t(MposWithEdges withEdges = MposWithEdges::OFF);
    std::vector<Eigen::Tensor<cplx, 4>>                get_compressed_mpos(MposWithEdges withEdges = MposWithEdges::OFF);
    std::vector<Eigen::Tensor<cplx, 4>>                get_compressed_mpos_squared(MposWithEdges withEdges = MposWithEdges::OFF);

    [[nodiscard]] bool                  has_energy_shifted_mpo() const; // For shifted energy MPO's
    [[nodiscard]] bool                  has_energy_shifted_mpo_squared() const; // For shifted energy MPO's
    [[nodiscard]] bool                  is_compressed_mpo_squared() const;
    [[nodiscard]] double                get_energy_shift() const;
    [[nodiscard]] double                get_energy_shift_per_site() const;
    [[nodiscard]] std::vector<std::any> get_parameter(std::string_view fieldname);

    // For local operations
    ModelLocal get_local(const std::vector<size_t> &sites) const;
    ModelLocal get_local() const;

    // For multisite
    std::array<long, 4>    active_dimensions() const;
    Eigen::Tensor<cplx, 4> get_multisite_mpo(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt, bool with_edgeL = false,
                                             bool with_edgeR = false) const;
    Eigen::Tensor<cplx_t, 4>        get_multisite_mpo_t(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt,
                                                        bool with_edgeL = false, bool with_edgeR = false) const;
    Eigen::Tensor<cplx, 2>          get_multisite_ham(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    Eigen::Tensor<cplx_t, 2>        get_multisite_ham_t(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    const Eigen::Tensor<cplx, 4>   &get_multisite_mpo() const;
    const Eigen::Tensor<cplx, 2>   &get_multisite_ham() const;
    const Eigen::Tensor<cplx_t, 4> &get_multisite_mpo_t() const;
    const Eigen::Tensor<cplx_t, 2> &get_multisite_ham_t() const;

    Eigen::Tensor<cplx, 4> get_multisite_mpo_shifted_view(double energy_per_site) const;
    Eigen::Tensor<cplx, 4> get_multisite_mpo_squared_shifted_view(double energy_per_site) const;

    std::array<long, 4>           active_dimensions_squared() const;
    Eigen::Tensor<cplx, 4>        get_multisite_mpo_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    Eigen::Tensor<cplx, 2>        get_multisite_ham_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody = std::nullopt) const;
    const Eigen::Tensor<cplx, 4> &get_multisite_mpo_squared() const;
    const Eigen::Tensor<cplx, 2> &get_multisite_ham_squared() const;
    void                          clear_cache(LogPolicy logPolicy = LogPolicy::QUIET) const;

    std::vector<size_t> get_active_ids() const;
    std::vector<size_t> get_active_ids_sq() const;
};
