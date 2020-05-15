//
// Created by david on 2020-05-12.
//

#include "class_env_base.h"
#include <tensors/model/class_mpo_base.h>
#include <tensors/state/class_mps_site.h>

using Scalar = class_env_base::Scalar;

class_env_base::class_env_base(std::string side_, size_t position_):position(position_),side(side_){}
class_env_base::class_env_base(std::string side_, const class_mps_site & MPS, const class_mpo_base &MPO):side(side_)
{
    if (MPS.get_position() != MPO.get_position())
        throw std::logic_error(fmt::format("MPS and MPO have different positions: {} != {}", MPS.get_position(),MPO.get_position()));
    position = MPS.get_position();
}
