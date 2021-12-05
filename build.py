import os
import errno
import warnings
import argparse
import subprocess
from packaging import version
import re
import glob
import shutil
import psutil

def str2bool(v):
    if isinstance(v, bool):
            return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse(project_name):
    parser = argparse.ArgumentParser(description='CMake Project Builder for {}'.format(project_name))
    parser.add_argument('-a', '--arch', type=str, help='Choose microarchitecture',
                        default='haswell', choices=['core2', 'nehalem', 'sandybridge', 'haswell', 'native'])
    parser.add_argument('-b', '--build-type', type=str, help='Build type',
                        default='Release', choices=['Release', 'Debug', 'RelWithDebInfo'])
    parser.add_argument('--clear-cmake', action='store_true', help='Delete CMake build directory before build')
    parser.add_argument('-c', '--clear-cache', action='store_true', help='Delete CMakeCache.txt before build')
    parser.add_argument('-d', '--dry-run', action='store_true', help='Dry run makes no changes')
    parser.add_argument('-e', '--examples', action='store_true', help='Build examples')
    parser.add_argument('-i', '--install-prefix', type=str, help='Install Prefix', default='install')
    parser.add_argument('-j', '--make-threads', type=int, help='Make Threads', default=psutil.cpu_count(logical = False))
    parser.add_argument('-G', '--generator', type=str, help='CMake Generator', default=None,
                        choices=[None, 'Ninja', 'Unix Makefiles'])
    parser.add_argument('-p', '--package-manager', type=str, help='Package Manager', default='conan',
                        choices=['find', 'cmake', 'find-or-cmake', 'conan'])
    parser.add_argument('-s', '--shared', type=str2bool, nargs='?', const=True, default=False, help="Shared library linkage")
    parser.add_argument('-t', '--target', type=str, help='Build Target', default=None)
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose build')
    parser.add_argument('-X', '--cmake-args', action='append', type=str, help='Extra CMake arguments', default=[])
    parser.add_argument('--mkl', action='store_true', help='Use Intel Math Kernel Library as instead of OpenBLAS')
    parser.add_argument('--coverage', action='store_true', help='Enable test Coverage')
    parser.add_argument('--lto', action='store_true', help='Enable Link Time Optimization')
    parser.add_argument('--pch', action='store_true', help='Enable Precompiled Headers')
    parser.add_argument('--ccache', action='store_true', help='Enable ccache to speed up repeated builds')
    parser.add_argument('--asan', action='store_true', help='Enable Address Sanitizer')
    parser.add_argument('--test', action='store_true', help='Enable CTest. Build and run tests')
    parser.add_argument('--threads', action='store_true', help='Enable STL Threading in Eigen::Tensor')
    parser.add_argument('--loglevel', type=str, help='CMake Log Level', default=None,
                        choices=['TRACE', 'DEBUG', 'VERBOSE', 'STATUS', 'NOTICE'])
    parser.add_argument('--build-dir', type=str, help='CMake build directory', default='build')
    parser.add_argument('--default-kraken', action='store_true', help='Set defaults for kraken cluster')
    parser.add_argument('--default-tetralith', action='store_true', help='Set defaults for tetralith cluster')
    parser.add_argument('--default-desktop', action='store_true', help='Set defaults for a regular desktop')
    parser.add_argument('--default-actions', action='store_true', help='Set defaults for a GitHub Actions')
    parser.add_argument('--disable-color', action='store_true', help='Disable color output')
    parser.add_argument('--debug', action='store_true', help='Debug this builder script')


    args = parser.parse_args()
    if args.default_kraken:
        parser.set_defaults(mkl=True, lto=True, pch=True, ccache=True, generator='Ninja', package_manager='conan', arch='haswell')
        args = parser.parse_args()
    if args.default_tetralith:
        parser.set_defaults(mkl=True, lto=True, pch=True, ccache=True, generator='Ninja', package_manager='conan', arch='native', make_threads=16)
        args = parser.parse_args()
    if args.default_desktop:
        parser.set_defaults(mkl=True, lto=True, pch=True, ccache=True, generator='Ninja', package_manager='conan', arch='native')
        args = parser.parse_args()
    if args.default_desktop:
        parser.set_defaults(mkl=False, lto=True, pch=True, ccache=True, generator=None, package_manager='conan', arch='generic', verbose=True, make_threads=2, test=True, asan=True, coverage=True, )
        args = parser.parse_args()
    return args


def assert_cmake_version(cmake_min_version):
    cmake_which_cmd = subprocess.Popen(['which', 'cmake'], shell=False, stdout=subprocess.PIPE, encoding='utf-8')
    which,_ = cmake_which_cmd.communicate()
    cmake_which_cmd.stdout.close()
    cmake_version_cmd = subprocess.Popen(['cmake', '--version'], shell=False, stdout=subprocess.PIPE, encoding='utf-8')
    out, err = cmake_version_cmd.communicate()
    cmake_version_cmd.stdout.close()

    if not err:
        cmake_version = version.parse(re.search(r'version \s*([\d.]+)', out).group(1))
        if cmake_version < cmake_min_version:
            raise ValueError(
                'CMake version {} is lower than the minimum required: {}'.format(cmake_version, cmake_min_version))
        print('Found CMake version {} at {}'.format(cmake_version, which))

def silentrm(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

def silentrmtree(path):
    try:
        shutil.rmtree(path)
    except:
        print("Directory could not be removed:", path)
        pass


def generate_cmake_commands(project_name, args):
    cmake_cfg = ['cmake']
    cmake_bld = ['cmake']
    cmake_tst = []
    cmake_env = os.environ.copy()  # Add environment variables here
    cmake_c_flags_init = []
    cmake_cxx_flags_init = []


    cmake_cfg.extend(['-S', '.'])  # Path to source
    cmake_cfg.extend(['-B', '{}/{}'.format(args.build_dir, args.build_type)])  # Path to build
    cmake_bld.extend(['--build', '{}/{}'.format(args.build_dir, args.build_type)])


    if not args.disable_color:
        cmake_env['CLICOLOR'] = str(1)
        cmake_env['CLICOLOR_FORCE'] = str(1)
        cmake_c_flags_init.extend(['-fdiagnostics-color=always'])
        cmake_cxx_flags_init.extend(['-fdiagnostics-color=always'])

    if args.arch:
        cmake_cfg.extend(['-D{}_MICROARCH:STRING={}'.format(project_name.upper(), args.arch)])

    if args.build_type:
        cmake_cfg.extend(['-DCMAKE_BUILD_TYPE:STRING={}'.format(args.build_type)])

    if args.clear_cache:
        print('Deleting {}/{}/CMakeCache.txt'.format(args.build_dir, args.build_type))
        silentrm('{}/{}/CMakeCache.txt'.format(args.build_dir, args.build_type))
        conanFiles = glob.glob('{}/{}/conan*'.format(args.build_dir, args.build_type))
        for file in conanFiles:
            silentrm(file)

    if args.clear_cmake:
        print('Deleting directory {}/{}'.format(args.build_dir, args.build_type))
        silentrmtree('{}/{}'.format(args.build_dir, args.build_type))

    if args.install_prefix:
        cmake_cfg.extend(['-DCMAKE_INSTALL_PREFIX:PATH={}'.format(args.install_prefix)])

    if args.make_threads:
        cmake_env['CMAKE_BUILD_PARALLEL_LEVEL'] = str(args.make_threads)
        cmake_bld.extend(['--parallel', str(args.make_threads)])

    if args.generator:
        cmake_env['CMAKE_GENERATOR'] = args.generator
        cmake_cfg.extend(['-G', args.generator])

    if args.package_manager:
        cmake_cfg.extend(['-D{}_PACKAGE_MANAGER:STRING={}'.format(project_name.upper(), args.package_manager)])

    if args.shared:
        cmake_cfg.extend(['-DBUILD_SHARED_LIBS:BOOL=ON'])

    if args.target:
        cmake_bld.extend(['--target', args.target])

    if args.verbose:
        cmake_cfg.extend(['-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON'])

    if args.cmake_args:
        cmake_cfg.extend(args.cmake_args)

    if args.mkl:
        cmake_cfg.extend(['-D{}_ENABLE_MKL:BOOL=ON'.format(project_name.upper())])

    if args.coverage:
        cmake_cfg.extend(['-D{}_ENABLE_COVERAGE:BOOL=ON'.format(project_name.upper())])

    if args.lto:
        cmake_cfg.extend(['-D{}_ENABLE_LTO:BOOL=ON'.format(project_name.upper())])

    if args.pch:
        cmake_cfg.extend(['-D{}_ENABLE_PCH:BOOL=ON'.format(project_name.upper())])

    if args.ccache:
        cmake_cfg.extend(['-D{}_ENABLE_CCACHE:BOOL=ON'.format(project_name.upper())])

    if args.asan:
        cmake_cfg.extend(['-D{}_ENABLE_ASAN:BOOL=ON'.format(project_name.upper())])

    if cmake_c_flags_init:
        cmake_cfg.extend(['-DCMAKE_C_FLAGS_INIT:STRING={}'.format(' '.join(cmake_c_flags_init))])
    if cmake_cxx_flags_init:
        cmake_cfg.extend(['-DCMAKE_CXX_FLAGS_INIT:STRING={}'.format(' '.join(cmake_cxx_flags_init))])

    if args.test:
        cmake_cfg.extend(['-D{}_ENABLE_TESTS:BOOL=ON'.format(project_name.upper())])
        cmake_tst.extend(['ctest', '--output-on-failure', '--build-config', args.build_type,
                          '--test-dir', '{}/{}'.format(args.build_dir, args.build_type)])
        if args.target and 'test' in args.target:
            cmake_tst.extend(['--build-target', args.target, '-R', args.target])  # Pattern match only the built target
        else:
            cmake_tst.extend(['-R', project_name.lower()])  # Pattern match test targets containing the project name

    if args.threads:
        cmake_cfg.extend(['-D{}_ENABLE_THREADS:BOOL=ON'.format(project_name.upper())])

    if args.examples:
        cmake_cfg.extend(['-D{}_BUILD_EXAMPLES:BOOL=ON'.format(project_name.upper())])

    if args.loglevel:
        cmake_cfg.extend(['--loglevel={}'.format(args.loglevel)])

    if args.debug:
        print(' '.join(cmake_cfg))
        print(' '.join(cmake_bld))
        print(' '.join(cmake_tst))
        print(' '.join(cmake_env))
    return cmake_cfg, cmake_bld, cmake_tst, cmake_env

def run(cmd,env, args):
    if len(cmd) == 0:
        return
    if args.debug:
        print("Running command: ", ' '.join(cmd))
    with subprocess.Popen(cmd,bufsize=1, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8',env=env) as p:
        for line in p.stdout:
            print(line, end='')  # process line here
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, ' '.join(p.args))


def main():
    project_name = 'DMRG'
    cmake_min_version = version.parse("3.18")
    assert_cmake_version(cmake_min_version=cmake_min_version)
    args = parse(project_name=project_name)
    cmake_cfg_cmd, cmake_bld_cmd, cmake_tst_cmd, cmake_env = generate_cmake_commands(project_name=project_name,
                                                                                     args=args)
    print("Building with {} threads".format(args.make_threads))
    run(cmd=cmake_cfg_cmd, env=cmake_env, args=args)
    run(cmd=cmake_bld_cmd, env=cmake_env, args=args)
    run(cmd=cmake_tst_cmd, env=cmake_env, args=args)



if __name__ == "__main__":
    main()
