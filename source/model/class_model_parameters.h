//
// Created by david on 2018-10-01.
//

#pragma once

#include <array>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>

struct p_selfdual_tf_rf_ising {
    double J_rnd            = 0;         /*!< Randomly distributed nearest neighbour coupling */
    double h_rnd            = 0;         /*!< Randomly distributed on-site field */
    double J_ptb            = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
    double h_ptb            = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
    double J_avg            = 0;         /*!< Average of J_rnd between all sites*/
    double h_avg            = 0;         /*!< Average of h_rnd on all sites */
    double J_mean           = 0;         /*!< Mean for the distrbution of J_rnd */
    double h_mean           = 0;         /*!< Mean for the distrbution of h_rnd */
    double J_sigma          = 0;         /*!< Standard deviation for the distribution of J_rnd */
    double h_sigma          = 0;         /*!< Standard deviation for the distribution of h_rnd */
    double lambda           = 0;         /*!< Factor involved in next-nearest neighbor interaction */
    double delta            = 0;         /*!< Difference J_mean - h_mean  */
    size_t spin_dim         = 0;         /*!< Spin dimension */
    char   distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
};

struct p_tf_ising {
    double J_nn             = 0;         /*!< Nearest neighbor coupling */
    double J_nnn            = 0;         /*!< Next-nearest neighbor coupling */
    double h_field          = 0;         /*!< On-site magnetic field */
    double h_rnd            = 0;         /*!< Random field value */
    double h_ptb            = 0;         /*!< Perturbation */
    double h_mean           = 0;         /*!< Mean of the distribution for the random field */
    double h_sigma          = 0;         /*!< Randomness strength. In distribution this is N(r_mean,r_sigma) or U(r_mean-r_sigma,r_mean+r_sigma) */
    size_t spin_dim         = 0;         /*!< Spin dimension */
    char   distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
};
