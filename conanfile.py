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
        if self.settings.compiler != "gcc" and self.options.with_quadmath:
            self.output.warning("Quadmath disabled (requires a GCC compiler)")
            self.options.with_quadmath = False

        self.options["hdf5"].with_zlib = self.options.with_zlib
        self.options["h5pp"].with_quadmath = self.options.with_quadmath
        self.options["ceres-solver"].use_glog = True
        self.options["ceres-solver"].use_CXX11_threads = True
        self.options["ceres-solver"].use_eigen_sparse = False

        # Use backtrace because this doesn't invoke any further dependencies that break easily on CI.
        # We take care of linking with system libdw instead, if needed.
        self.options["backward-cpp"].stack_walking = "backtrace"
        self.options["backward-cpp"].stack_details = "backtrace_symbol"



    def requirements(self):
        self.requires("h5pp/[>=1.11.3 <1.12]@davidace/dev")
        self.requires("fmt/[>=10.2.0 <11]")
        self.requires("spdlog/[>=1.14.0 <1.15]")
        self.requires("eigen/3.4.0")
        self.requires("ceres-solver/2.2.0")
        self.requires("cli11/2.3.2")
        self.requires("backward-cpp/1.6")
        self.requires("pcg-cpp/cci.20220409")
        self.requires("tomlplusplus/3.4.0")



    def validate(self):
        if self.settings.compiler.get_safe("cppstd"):
            check_min_cppstd(self, 17)
        minimum_version = self._compilers_minimum_version.get(str(self.settings.compiler), False)
        if minimum_version:
            if Version(self.settings.compiler.version) < minimum_version:
                raise ConanInvalidConfiguration("DMRG++ requires C++17, which your compiler does not support.")
        else:
            self.output.warning("DMRG++ requires C++17. Your compiler is unknown. Assuming it supports C++17.")


### NOTE October 4, 2020 ####
# Another flag that seems to fix weird release-only bugs is
#            -fno-strict-aliasing

