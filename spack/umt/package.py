# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import socket
import itertools
from os import environ as env

from spack.package import *

class Umt(CachedCMakePackage, CudaPackage):
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

    variant("openmp", default=False, description="Enable OpenMP support")
    variant("openmp_offload", default=False, description="Enable OpenMP target offload support")
    variant("caliper", default=False, description="Enable Caliper performance timers")
    variant("umpire", default=False, description="Enable use of Umpire memory library")
    variant("shared", default=False, description="Enable shared libraries")
    variant("silo", default=False, description="Enable silo I/O support")
    variant("find_mpi", default=True, description="Use CMake find_package(mpi) logic.  Disable to rely on mpicxx, mpif90 compiler wrappers")

    conflicts('cuda_arch=none', when='+cuda', msg='CUDA architecture is required')
    ###########################################################################
    # package dependencies
    ###########################################################################

    #######################
    # CMake
    #######################
    depends_on("cmake@3.13.3:", type="build")

    #######################
    # Dependencies
    #######################
    depends_on("mpi", when="+find_mpi")
    depends_on("mpi+wrappers", when="~find_mpi")

    depends_on("cuda", when="+cuda")

    depends_on("conduit+fortran+shared", when="+shared")
    depends_on("conduit+fortran~shared", when="~shared")

    depends_on("mfem+conduit+mpi+shared", when="+mfem+shared")
    depends_on("mfem+conduit+mpi~shared", when="+mfem~shared")

    depends_on("hypre", when="+mfem")
    depends_on("metis", when="+mfem")

    depends_on("mfem+conduit+mpi~shared", when="+mfem~shared")
    depends_on("caliper+fortran+shared", when="+caliper+shared")
    depends_on("caliper+fortran~shared", when="+caliper~shared")

    depends_on("adiak", when="+caliper+shared")
    depends_on("adiak~shared", when="+caliper~shared")

    depends_on("umpire+fortran+shared", when="+umpire+shared")
    depends_on("umpire+fortran~shared", when="+umpire~shared")

    depends_on("camp", when="+umpire")

    depends_on("silo+shared", when="+silo+shared")
    depends_on("silo~shared", when="+silo~shared")

    ####################################################################
    # Note: cmake, build, and install stages are handled by CMakePackage
    ####################################################################

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

    # Override the parent class function until its fixed in spack.  (Greg Becker and Chris White are aware of bug)
    def flag_handler(self, name, flags):
        if name in ("cflags", "cxxflags", "cppflags", "fflags", "ldflags"):
            return (None, None, None)  # handled in the cmake cache
        return (flags, None, None)

    def initconfig_compiler_entries(self):
        spec = self.spec

        if "umpire" in spec:
            if "+rocm" in spec["umpire"]:
                spec.compiler_flags["cxxflags"].append("-I{0}/include".format(spec["hip"].prefix))

        #######################
        # Note - call the super class AFTER changing any flags, as the super class
        # adds the cflags, cxxflags, fflags, ldflags, etc, to the cache entries list.
        # If you try adding any of these yourself you will end up with duplicates.
        # - black27
        entries = super(Umt, self).initconfig_compiler_entries()

        if spec.satisfies("%cce"):
            entries.append(cmake_cache_option("STRICT_FPP_MODE", True))
            if "+openmp" in spec:
                entries.append(cmake_cache_option("OPENMP_HAS_USE_DEVICE_ADDR", True))
                entries.append(cmake_cache_option("OPENMP_HAS_FORTRAN_INTERFACE", True))

        if (len(self.compiler.extra_rpaths) > 0):
            # Provide extra link options to embed rpaths to libraries.
            # Spack is providing both the linker pass-through flag and the rpath flag in the
            # cc_rpath_arg string.  UMT CMake logic uses the target_link_options() command
            # to add these to its link and that requires just the paths.  Strip out the
            # linker pass through flags before handing to CMake.
            rpath_arg = self.compiler.cc_rpath_arg.replace(self.compiler.linker_arg, "")

            link_options = []
            link_options.extend( [rpath_arg + path for path in self.compiler.extra_rpaths] )
            entries.append(cmake_cache_string("TETON_LINK_OPTIONS", ",".join(link_options) ))

        return entries

    def initconfig_hardware_entries(self):
        spec = self.spec
        entries = super(Umt, self).initconfig_hardware_entries()

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

        return entries

    def initconfig_mpi_entries(self):
        entries = super(Umt, self).initconfig_mpi_entries()
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
        entries.append(cmake_cache_option("ENABLE_TESTS", True))

        found_hdf5_dependency = False
        found_zlib_dependency = False

        if "+silo" in self.spec:
            entries.append(cmake_cache_option("ENABLE_SILO", True))
            entries.append(cmake_cache_path("SILO_ROOT", self.spec["silo"].prefix))
            if "+hdf5" in spec["silo"]:
                found_hdf5_dependency = True
        else:
            entries.append(cmake_cache_option("ENABLE_SILO", False))

        entries.append(cmake_cache_path("CONDUIT_ROOT", spec["conduit"].prefix))

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

        if "+hdf5" in spec["conduit"]:
            found_hdf5_dependency = True

        if found_hdf5_dependency:
            entries.append(cmake_cache_option("ENABLE_HDF5", True))
            entries.append(cmake_cache_path("HDF5_ROOT", spec["hdf5"].prefix))

# Currently, UMT assumes that HDF5 was built with ZLIB.  Need to fix this.
            found_zlib_dependency = True
#            if "+zlib" in spec["hdf5"]:
#                found_zlib_dependency = True

        if "+zlib" in spec["conduit"]:
            found_zlib_dependency = True
        if found_zlib_dependency:
            entries.append(cmake_cache_path("Z_ROOT", spec["zlib"].prefix))

        return entries
