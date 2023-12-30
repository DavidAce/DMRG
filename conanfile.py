from conan import ConanFile
from conan.tools.layout import basic_layout
from conan.tools.files import copy
from conan.tools.microsoft import is_msvc
from conan.tools.build import check_min_cppstd
from conan.tools.scm import Version
from conan.errors import ConanInvalidConfiguration
import os

required_conan_version = ">=1.54.0"

class DMRGConan(ConanFile):
    name = "DMRG++"
    version = "3.1.0"
    description = "MPS algorithms for 1D quantum spin chains"
    homepage = "https://github.com/DavidAce/DMRG"
    author = "DavidAce <aceituno@kth.se>"
    topics = ("MPS", "DMRG", "TEBD")
    url = "https://github.com/DavidAce/DMRG"
    license = "MIT"
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps"
    no_copy_source = True
    short_paths = True
    exports_sources = "include/*"
    options = {
        "with_zlib" : [True, False],
        "with_quadmath": [True, False]
    }
    default_options = {
        "with_zlib" : True,
        "with_quadmath": True
    }

    @property
    def _compilers_minimum_version(self):
        return {
            "gcc": "7.4",
            "Visual Studio": "15.7",
            "clang": "6",
            "apple-clang": "10",
        }

    def config_options(self):
        self.options["hdf5"].with_zlib = self.options.with_zlib
        self.options["h5pp"].with_quadmath = self.options.with_quadmath
        self.options["ceres-solver"].use_glog = False
        # Use backtrace because this doesn't invoke any further dependencies that break easily on CI.
        # We take care of linking with system libdw instead, if needed.
        self.options["backward-cpp"].stack_walking = "backtrace"
        self.options["backward-cpp"].stack_details = "backtrace_symbol"
        if self.settings.compiler != "gcc":
            self.output.warning("Quadmath disabled (requires a GCC compiler)")
            self.options["h5pp"].with_quadmath=False


    def requirements(self):
        self.requires("h5pp/1.11.1@davidace/dev")
        self.requires("ceres-solver/2.2.0")
        self.requires("fmt/10.1.1")
        self.requires("spdlog/1.12.0")
        self.requires("eigen/3.4.0")
        self.requires("arpack++/2.3.0@davidace/dev")
        self.requires("cli11/2.3.2")
        self.requires("backward-cpp/1.6")
        self.requires("zlib/1.3")
        self.requires("pcg-cpp/cci.20220409")
        # self.requires("mpfr/4.1.0")


    def validate(self):
        if self.settings.compiler.get_safe("cppstd"):
            check_min_cppstd(self, 17)
        minimum_version = self._compilers_minimum_version.get(str(self.settings.compiler), False)
        if minimum_version:
            if Version(self.settings.compiler.version) < minimum_version:
                raise ConanInvalidConfiguration("DMRG++ requires C++17, which your compiler does not support.")
        else:
            self.output.warning("DMRG++ requires C++17. Your compiler is unknown. Assuming it supports C++17.")


### NOTE October 4 2020 ####
# Another flag that seems to fix weird release-only bugs is
#            -fno-strict-aliasing

### NOTE August 26 2020 ####
#
# Ceres started crashing on Tetralith again using -march=native.
# Tried to solve this issue once and for all.
# I've tried the following flags during compilation of DMRG++ and ceres-solver:
#
#           -DEIGEN_MALLOC_ALREADY_ALIGNED=[none,0,1]
#           -DEIGEN_MAX_ALIGN_BYTES=[none,16,32]
#           -march=[none,native]
#           -std=[none,c++17]
#
# Up until now, [0,16,none,none] has worked but now for some reason it stopped working.
# I noticed the -std=c++17 flag was not being passed on conan builds, so ceres defaulted to -std=c++14 instead.
# I fixed this in the conanfile.py of the ceres build. The -download-method=fetch method already had this fixed.
# When no Eigen flags were passed, and ceres-solver finally built with -std=c++17 the issues vanished.
# In the end what worked was [none,none,native,c++17] in both DMRG++ and ceres-solver.
# It is important that the same eigen setup is used in all compilation units, and c++17/c++14 seems to
# make Eigen infer some flags differently. In any case, settinc c++17 and no flags for eigen anywhere
# lets Eigen do its thing in the same way everywhere.
#
#
# ceres-solver:eigen_malloc_already_aligned=1
# ceres-solver:eigen_max_align_bytes=32
