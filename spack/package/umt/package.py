# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import socket
import itertools
from os import environ as env
import llnl.util.filesystem as fs

from spack.package import *

class Umt(CachedCMakePackage, CudaPackage, ROCmPackage):
    """Umt is a LLNL mini-app based on the Teton thermal radiative transport library."""

    homepage = "https://github.com/LLNL/UMT"
    url = ""
    git = 'https://github.com/LLNL/UMT.git'

    version("master", branch="master", submodules=False)
    maintainers = ["aaroncblack"]

    # The CMakeLists.txt is in 'src' directory.
    root_cmakelists_dir = "src"

    ###########################################################################
    # package variants
    ###########################################################################

    variant("fpp", default=False, description="Use simpler macros compatible with stricter Fortran preprocessors")

    variant("openmp", default=False, description="Enable OpenMP support")
    variant("openmp_offload", default=False, description="Enable OpenMP target offload support", when="+openmp")

    variant("caliper", default=False, description="Enable Caliper performance timers")
    variant("umpire", default=False, description="Enable use of Umpire memory library")
    variant("shared", default=False, description="Enable shared libraries")
    variant("silo", default=False, description="Enable silo I/O support")
    variant("find_mpi", default=True, description="Use CMake find_package(mpi) logic.  Disable to rely on mpicxx, mpif90 compiler wrappers")
    variant("tests", default=True, description="Enable test driver.")
    variant("host_config_only", default=False, description="Installs only the cmake cache file, for use in debugging.")

    conflicts('cuda_arch=none', when='+cuda', msg='CUDA architecture is required')
    ###########################################################################
    # package dependencies
    ###########################################################################

    #######################
    # CMake
    #######################
    depends_on("cmake@3.21.1:", type="build")

    #######################
    # Dependencies
    #######################
    depends_on("mpi")

    depends_on("cuda", when="+cuda")

    depends_on("conduit+fortran")
    depends_on("conduit+shared", when="+shared")
    depends_on("conduit~shared", when="~shared")

    depends_on("caliper+fortran", when="+caliper")
    depends_on("caliper+shared", when="+caliper+shared")
    depends_on("caliper~shared", when="+caliper~shared")

    depends_on("umpire+fortran", when="+umpire")
    depends_on("umpire+shared", when="+umpire+shared")
    depends_on("umpire~shared", when="+umpire~shared")
    depends_on("umpire+rocm", when="+rocm")
    depends_on("umpire+cuda", when="+cuda")
    depends_on("umpire+openmp", when="+openmp")

    depends_on("silo", when="+silo")
    depends_on("silo+shared", when="+silo+shared")
    depends_on("silo~shared", when="+silo~shared")

    ####################################################################
    # Note: cmake, build, and install stages are handled by CMakePackage
    ####################################################################

    def cmake(self, pkg, spec):
        if "+host_config_only" not in self.spec:
            super().cmake(pkg, spec)

    def build(self, pkg, spec):
        if "+host_config_only" not in self.spec:
            super().build(pkg, spec)

    def install(self, pkg, spec):
        if "+host_config_only" in self.spec:
            print ("Installing host config file only")
            fs.mkdirp(self.spec.prefix.share.cmake)
            fs.install(self.cache_path, self.spec.prefix.share.cmake)
        else:
            super().install(pkg, spec)

    def _get_sys_type(self, spec):
        sys_type = spec.architecture
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        return sys_type

    @property
    def cache_name(self):
        hostname = socket.gethostname()
        if "SYS_TYPE" in env:
            hostname = hostname.rstrip("1234567890")
        return "{0}-{1}-{2}@{3}.cmake".format(
            hostname,
            self._get_sys_type(self.spec),
            self.spec.compiler.name,
            self.spec.compiler.version,
        )

    def initconfig_compiler_entries(self):
        spec = self.spec

        # Spack is providing both the linker pass-through flag and the rpath flag in the
        # cc_rpath_arg string.  UMT CMake logic uses the target_link_options() command
        # to add these to its link and that requires just the paths.  Strip out the
        # linker pass through flags before handing to CMake.
        rpath_arg = self.compiler.cc_rpath_arg.replace(self.compiler.linker_arg, "")
        link_options = []

        #######################
        # Note - call the super class AFTER changing any flags, as the super class
        # adds the cflags, cxxflags, fflags, ldflags, etc, to the cache entries list.
        # If you try adding any of these yourself you will end up with duplicates.
        # - Aaron Black
        entries = super().initconfig_compiler_entries()

        if "+fpp" in spec:
            entries.append(cmake_cache_option("STRICT_FPP_MODE", True))

        if (len(self.compiler.extra_rpaths) > 0):
            # Provide extra link options to embed rpaths to libraries.

            link_options.extend( [rpath_arg + path for path in self.compiler.extra_rpaths] )
            entries.append(cmake_cache_string("TETON_LINK_OPTIONS", ",".join(link_options) ))

        return entries

    def initconfig_hardware_entries(self):
        spec = self.spec
        entries = super().initconfig_hardware_entries()

        #######################
        # Parallelism
        #######################
        if "+openmp" in spec:
            entries.append(cmake_cache_option("ENABLE_OPENMP", True))
        if "+openmp_offload" in spec:
            entries.append(cmake_cache_option("ENABLE_OPENMP_OFFLOAD", True))

        if "+cuda" in spec:
            entries.append(cmake_cache_option("ENABLE_CUDA", True))
            cuda_arch = spec.variants["cuda_arch"].value
            entries.append(cmake_cache_string("CMAKE_CUDA_ARCHITECTURES", "{0}".format(cuda_arch[0])))
            # Add CUDAToolkit_ROOT, as Spack does not set this.
            entries.append(cmake_cache_string("CUDAToolkit_ROOT", "{0}".format( spec["cuda"].prefix)))

        else:
            entries.append(cmake_cache_option("ENABLE_CUDA", False))

        if "+rocm" in spec:
            entries.append(cmake_cache_option("ENABLE_HIP", True))

        else:
            entries.append(cmake_cache_option("ENABLE_HIP", False))
        return entries

    def initconfig_mpi_entries(self):
        entries = super().initconfig_mpi_entries()
        if "+find_mpi" in self.spec:
            entries.append(cmake_cache_option("ENABLE_FIND_MPI", True))
        elif "~find_mpi" in self.spec:
            entries.append(cmake_cache_option("ENABLE_FIND_MPI", False))

        return entries

    def initconfig_package_entries(self):
        spec = self.spec
        entries = []

        #######################
        # Disable features not needed by UMT
        #######################
        entries.append(cmake_cache_option("ENABLE_MINIAPP_BUILD", True))


        if "+tests" in self.spec:
            entries.append(cmake_cache_option("ENABLE_TESTS", True))

        if "+silo" in self.spec:
            entries.append(cmake_cache_option("ENABLE_SILO", True))
            entries.append(cmake_cache_path("SILO_ROOT", self.spec["silo"].prefix))
        else:
            entries.append(cmake_cache_option("ENABLE_SILO", False))

        entries.append(cmake_cache_path("CONDUIT_ROOT", spec["conduit"].prefix))
        if "+parmetis" in spec["conduit"]:
            entries.append(cmake_cache_path("METIS_ROOT", spec["metis"].prefix))
            entries.append(cmake_cache_path("PARMETIS_ROOT", spec["parmetis"].prefix))
        if "+hdf5" in spec["conduit"]:
            need_hdf5 = True
        if "+zlib" in spec["conduit"]:
            need_zlib = True

        if "+caliper" in spec:
            entries.append(cmake_cache_option("ENABLE_CALIPER", True))
            entries.append(cmake_cache_path("CALIPER_ROOT", spec["caliper"].prefix))
            if "adiak" in spec:
                entries.append(cmake_cache_path("ADIAK_ROOT", spec["adiak"].prefix))

        if "+umpire" in spec:
            entries.append(cmake_cache_option("ENABLE_UMPIRE", True))
            entries.append(cmake_cache_path("UMPIRE_ROOT", spec["umpire"].prefix))
            entries.append(cmake_cache_option("ENABLE_CAMP", True))
            entries.append(cmake_cache_path("CAMP_ROOT", spec["camp"].prefix))
            entries.append(cmake_cache_option("ENABLE_FMT", True))
            entries.append(cmake_cache_path("FMT_ROOT", spec["fmt"].prefix))

        # Silo or Conduit may pull in HDF5
        if "+hdf5" in spec["conduit"] or ("+silo" in self.spec and "+hdf5" in spec["silo"]):
            entries.append(cmake_cache_option("ENABLE_HDF5", True))
            entries.append(cmake_cache_path("HDF5_ROOT", spec["hdf5"].prefix))

            # HDF5 in turn depends on zlib
            if "zlib-api" in spec:
                entries.append(cmake_cache_path("Z_ROOT", spec["zlib-api"].prefix))

        return entries
