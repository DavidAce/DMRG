function(check_ceres_compiles TAG REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES REQUIRED_TARGETS)
    if(REQUIRED_LIBRARIES_UNPARSED AND REQUIRED_TARGETS)
        message(FATAL_ERROR "Please use EITHER required libs or required targets to honor linking order")
    endif()

    if(CMAKE_VERBOSE_MAKEFILE)
        message(STATUS "Checking if ceres compiles")
        message(STATUS "CERES TEST [${TAG}] REQUIRED_FLAGS                 ${REQUIRED_FLAGS}")
        message(STATUS "CERES TEST [${TAG}] REQUIRED_DEFINITIONS           ${REQUIRED_DEFINITIONS}")
        message(STATUS "CERES TEST [${TAG}] REQUIRED_LIBRARIES_UNPARSED    ${REQUIRED_LIBRARIES_UNPARSED}")
        message(STATUS "CERES TEST [${TAG}] REQUIRED_INCLUDES              ${REQUIRED_INCLUDES}")
        message(STATUS "CERES TEST [${TAG}] REQUIRED_TARGETS               ${REQUIRED_TARGETS}")
    endif()


    if(NOT BUILD_SHARED_LIBS)
        list(APPEND REQUIRED_LIBRARIES -static)
    endif()

    include(cmake-modules/getExpandedTarget.cmake)
    foreach(elem ${REQUIRED_LIBRARIES_UNPARSED})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            expand_target_incs(${elem} expanded_incs)
            expand_target_opts(${elem} expanded_opts)
            list(APPEND REQUIRED_LIBRARIES "${expanded_libs}")
            list(APPEND REQUIRED_INCLUDES  "${expanded_incs}")
            list(APPEND REQUIRED_FLAGS     "${expanded_opts}")
        else()
            list(APPEND REQUIRED_LIBRARIES "${elem}")
        endif()
    endforeach()

    foreach(elem ${REQUIRED_TARGETS})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            expand_target_incs(${elem} expanded_incs)
            expand_target_opts(${elem} expanded_opts)
            list(APPEND REQUIRED_LIBRARIES "${expanded_libs}")
            list(APPEND REQUIRED_INCLUDES  "${expanded_incs}")
            list(APPEND REQUIRED_FLAGS     "${expanded_opts}")
        endif()
    endforeach()
    string(REPLACE ";" " " REQUIRED_FLAGS      "${REQUIRED_FLAGS}") # Needs to be a space-separated list
    list(FILTER REQUIRED_LIBRARIES EXCLUDE REGEX "NOTFOUND") # Make sure every item is valid
    list(FILTER REQUIRED_INCLUDES  EXCLUDE REGEX "NOTFOUND")

    #   Test features
    set(CMAKE_REQUIRED_FLAGS        ${REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_DEFINITIONS  ${REQUIRED_DEFINITIONS})
    set(CMAKE_REQUIRED_LIBRARIES    ${REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES     ${REQUIRED_INCLUDES})
    if(CMAKE_VERBOSE_MAKEFILE)
        message(STATUS "CERES TEST [${TAG}] CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "CERES TEST [${TAG}] CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
        message(STATUS "CERES TEST [${TAG}] CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "CERES TEST [${TAG}] CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")
    endif()

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
        " CERES_COMPILES_${TAG})
    if(CERES_COMPILES_${TAG})
        set(CERES_COMPILES_${TAG} TRUE PARENT_SCOPE)
    else()
        set(CERES_COMPILES_${TAG} FALSE PARENT_SCOPE)
    endif()
endfunction()
