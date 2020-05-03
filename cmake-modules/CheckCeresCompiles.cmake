
function(check_ceres_compiles TARGETS LIBS INCS OPTS DEFS)
    if(NOT BUILD_SHARED_LIBS)
        list(APPEND CMAKE_REQUIRED_LIBRARIES -static)
    endif()

    list(APPEND CMAKE_REQUIRED_LIBRARIES     ${LIBS} ${TARGETS})
    list(APPEND CMAKE_REQUIRED_INCLUDES      ${INCS})
    list(APPEND CMAKE_REQUIRED_FLAGS         ${OPTS} -std=c++17)
    list(APPEND CMAKE_REQUIRED_DEFINITIONS   ${DEFS})

    list(TRANSFORM "CMAKE_REQUIRED_DEFINITIONS" PREPEND "-D")  # Definitions should start with "-D"
    string(REPLACE ";" " "  CMAKE_REQUIRED_FLAGS          "${CMAKE_REQUIRED_FLAGS}")        # Needs to be a space-separated list

    if(DMRG_PRINT_CHECKS)
        message(STATUS "CERES COMPILE TEST CMAKE_REQUIRED_LIBRARIES         ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "CERES COMPILE TEST CMAKE_REQUIRED_INCLUDES          ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "CERES COMPILE TEST CMAKE_REQUIRED_FLAGS             ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "CERES COMPILE TEST CMAKE_REQUIRED_DEFINITIONS       ${CMAKE_REQUIRED_DEFINITIONS}")
    endif()

    #   Test features
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
            #include <ceres/ceres.h>
            #include <glog/logging.h>
            // f(x,y) = (1-x)^2 + 100(y - x^2)^2;
            class Rosenbrock : public ceres::FirstOrderFunction {
             public:
              virtual ~Rosenbrock() {}
              virtual bool Evaluate(const double* parameters,
                                    double* cost,
                                    double* gradient) const {
                const double x = parameters[0];
                const double y = parameters[1];
                cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
                if (gradient != nullptr) {
                  gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
                  gradient[1] = 200.0 * (y - x * x);
                }
                return true;
              }
              virtual int NumParameters() const { return 2; }
            };
            int main(int argc, char** argv) {
              google::InitGoogleLogging(argv[0]);
              double parameters[2] = {-1.2, 1.0};
              ceres::GradientProblemSolver::Options options;
              options.minimizer_progress_to_stdout = true;
              ceres::GradientProblemSolver::Summary summary;
              ceres::GradientProblem problem(new Rosenbrock());
              ceres::Solve(options, problem, parameters, &summary);
              std::cout << summary.FullReport() << std::endl;
              return 0;
            }
        " CERES_COMPILES)
    set(CERES_COMPILES ${CERES_COMPILES} PARENT_SCOPE)
    if(NOT CERES_COMPILES)
        unset(CERES_COMPILES CACHE)
        unset(CERES_COMPILES PARENT_SCOPE)
        if(DMRG_PRINT_CHECKS AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
            file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
            message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
        endif()
    endif()
endfunction()
