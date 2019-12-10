//
// Created by david on 2019-12-10.
//

#pragma once
#include <tools/finite/opt.h>
#include "ceres/first_order_function.h"

namespace tools::finite::opt::internal {
        class ceres_rosenbrock_functor : public ceres::FirstOrderFunction {
        public:
            bool Evaluate(const double *parameters,
                                  double *cost,
                                  double *gradient) const final
            {
                const double x = parameters[0];
                const double y = parameters[1];

                cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
                if (gradient != nullptr) {
                    gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
                    gradient[1] = 200.0 * (y - x * x);
                }
                return true;
            }

            int NumParameters() const final { return 2; }
        };
    }


