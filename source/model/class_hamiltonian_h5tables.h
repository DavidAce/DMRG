//
// Created by david on 2018-10-01.
//

#ifndef DMRG_SELFDUAL_TF_RF_ISING_H5log_H
#define DMRG_SELFDUAL_TF_RF_ISING_H5log_H

#include <vector>
#include <array>
#include <hdf5.h>
#include <hdf5_hl.h>




//Profiling table definition
class class_selfdual_tf_rf_ising_h5table{
private:
    struct data_struct{
        int    position;
        double J_rnd;
        double h_rnd;
        double J_log_mean;
        double h_log_mean;
        double J_sigma;
        double h_sigma;
        double lambda;
        double e_reduced;
        double spin_dim;


        data_struct(int    position_,
                    double J_rnd_     ,  double h_rnd_,
                    double J_log_mean_,  double h_log_mean_,
                    double J_sigma_   ,  double h_sigma_,
                    double lambda_    ,  double e_reduced_,
                    double spin_dim_)
                :position(position_),
                 J_rnd(J_rnd_)          , h_rnd(h_rnd_),
                 J_log_mean(J_log_mean_), h_log_mean(h_log_mean_),
                 J_sigma(J_sigma_)      , h_sigma(h_sigma_),
                 lambda(lambda_)        , e_reduced(e_reduced_),
                 spin_dim(spin_dim_)
        {}
    };
    struct meta_struct{
        constexpr static hsize_t                NFIELDS     = 10;
        size_t           dst_size                           = sizeof (data_struct);
        std::array       <size_t,NFIELDS>       dst_offsets = {HOFFSET(data_struct, position),
                                                               HOFFSET(data_struct, J_rnd     ), HOFFSET(data_struct,  h_rnd     ),
                                                               HOFFSET(data_struct, J_log_mean), HOFFSET(data_struct,  h_log_mean),
                                                               HOFFSET(data_struct, J_sigma   ), HOFFSET(data_struct,  h_sigma   ),
                                                               HOFFSET(data_struct, lambda    ), HOFFSET(data_struct,  e_reduced ),
                                                               HOFFSET(data_struct, spin_dim)
        };
        std::array       <size_t,NFIELDS>       dst_sizes   = {sizeof(data_struct::position),
                                                               sizeof(data_struct::J_rnd     ), sizeof(data_struct::h_rnd     ),
                                                               sizeof(data_struct::J_log_mean), sizeof(data_struct::h_log_mean),
                                                               sizeof(data_struct::J_sigma   ), sizeof(data_struct::h_sigma   ),
                                                               sizeof(data_struct::lambda    ), sizeof(data_struct::e_reduced ),
                                                               sizeof(data_struct::spin_dim)
        };
        std::array       <const char*,NFIELDS>  field_names = {"position",
                                                               "J_rnd"     , "h_rnd",
                                                               "J_log_mean", "h_log_mean",
                                                               "J_sigma"   , "h_sigma",
                                                               "lambda"    , "e_reduced",
                                                               "spin_dim"
        };

        std::array       <hid_t,NFIELDS>        field_types = {H5T_NATIVE_INT,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE
        };

        hsize_t          chunk_size                         = 1;
        void             *fill_data                         = nullptr;
        int              compress                           = 0;
    };
public:
    class_selfdual_tf_rf_ising_h5table() = default;
};



#endif //DMRG_log_SELFDUAL_TF_RF_ISING_H
