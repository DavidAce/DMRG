#pragma once
#include <array>
#include <complex>
#include <vector>

class MpoSite;
class TensorsFinite;
// class StateLocal;

class ModelLocal {
    public:
    using cplx = std::complex<double>;
    using real = double;

    std::vector<std::unique_ptr<MpoSite>> mpos; /*!< A subset of mpos */

    std::vector<size_t>              get_positions() const;
    std::vector<std::array<long, 4>> get_dimensions() const;

    //    operator * ()
};