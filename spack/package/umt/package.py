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

    version("develop", branch="develop", submodules=False)
    maintainers = ["aaroncblack"]

    # The CMakeLists.txt is in 'src' directory.
    root_cmakelists_dir = "src"

    ###########################################################################
    # package variants
    ###########################################################################

    variant("openmp", default=False, description="Enable OpenMP support")
    variant("openmp_offload", default=False, description="Enable OpenMP target offload support")
    variant("caliper", default=False, description="Enable Caliper performance timers")
    variant("umpire", default=False, description="Enable use of Umpire memory library")
    variant("find_mpi", default=True, description="Use CMake find_package(mpi) logic.  Disable to rely on mpicxx, mpif90 compiler wrappers")
    variant("tests", default=True, description="Enable test driver.")

    conflicts('cuda_arch=none', when='+cuda', msg='CUDA architecture is required')
    conflicts('amdgpu_target=none', when='+rocm', msg='AMD GPU architecture is required')

    ###########################################################################
    # package dependencies
    ###########################################################################
    depends_on("cmake@3.13.3:", type="build")

    depends_on("mpi", when="+find_mpi")
    depends_on("mpi+wrappers", when="~find_mpi")

    depends_on("cuda", when="+cuda")
    depends_on("hip", when="+rocm")

    depends_on("conduit+fortran")
    depends_on("caliper+fortran", when="+caliper")
    depends_on("umpire+fortran", when="+umpire")

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
        # - aaroncblack
        entries = super().initconfig_compiler_entries()

        if spec.satisfies("%cce"):
            entries.append(cmake_cache_option("STRICT_FPP_MODE", True))
            if "+openmp" in spec:
                entries.append(cmake_cache_option("OPENMP_HAS_USE_DEVICE_ADDR", True))
                entries.append(cmake_cache_option("OPENMP_HAS_FORTRAN_INTERFACE", True))

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
        else:
            entries.append(cmake_cache_option("ENABLE_FIND_MPI", False))

        return entries

    def initconfig_package_entries(self):
        spec = self.spec
        entries = []

        #######################
        # Disable features not needed by UMT
        #######################
        entries.append(cmake_cache_option("ENABLE_MINIAPP_BUILD", True))

        # Enable importing cmake targets.
        entries.append(cmake_cache_option("ENABLE_FIND_PACKAGE_CONFIG_MODE", True))

        if "+tests" in self.spec:
            entries.append(cmake_cache_option("ENABLE_TESTS", True))

        entries.append(cmake_cache_path("CONDUIT_ROOT", spec["conduit"].prefix))
        if "+parmetis" in spec["conduit"]:
            entries.append(cmake_cache_path("METIS_ROOT", spec["metis"].prefix))
            entries.append(cmake_cache_path("PARMETIS_ROOT", spec["parmetis"].prefix))
        if "+hdf5" in spec["conduit"]:
            need_hdf5 = True
        if "+zlib" in spec["conduit"]:
            entries.append(cmake_cache_path("Z_ROOT", spec["zlib"].prefix))

        if "+caliper" in spec:
            entries.append(cmake_cache_option("ENABLE_CALIPER", True))
            entries.append(cmake_cache_path("CALIPER_ROOT", spec["caliper"].prefix))
            if "+adiak" in spec["caliper"]:
                entries.append(cmake_cache_path("ADIAK_ROOT", spec["adiak"].prefix))

        if "+umpire" in spec:
            entries.append(cmake_cache_option("ENABLE_UMPIRE", True))
            entries.append(cmake_cache_path("UMPIRE_ROOT", spec["umpire"].prefix))
            entries.append(cmake_cache_option("ENABLE_CAMP", True))
            entries.append(cmake_cache_path("CAMP_ROOT", spec["camp"].prefix))
            if ("+fmt" in spec["umpire"]):
                entries.append(cmake_cache_option("ENABLE_FMT", True))
                entries.append(cmake_cache_path("FMT_ROOT", spec["fmt"].prefix))

        return entries
