# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import socket
from os import environ as env

from spack.package import *

class Umt(CachedCMakePackage):
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

    variant("openmp", default=True, description="Enable OpenMP support")
    variant("openmp_offload", default=False, description="Enable OpenMP target offload support")
    variant("caliper", default=False, description="Enable Caliper performance timers")
    variant("mpi", default=True, description="Enable MPI support (mandatory")
    variant("mfem", default=False, description="Enable support for reading MFEM meshes")
    variant("umpire", default=False, description="Enable use of Umpire memory library")

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
    depends_on("mpi")
    depends_on("conduit+fortran~hdf5")
    depends_on("mfem", when="+mfem")
    depends_on("caliper+fortran", when="+caliper")
    depends_on("adiak", when="+caliper")
    depends_on("umpire+fortran", when="+umpire")
    depends_on("camp", when="+umpire")

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

    def cmake_args(self):
        options = []
        return options

    def _get_host_config_path(self, spec):
        sys_type = spec.architecture
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        host_config_path = "{0}-{1}-{2}-{3}.cmake".format(
            socket.gethostname(), sys_type, spec.compiler, spec.dag_hash()
        )
        dest_dir = spec.prefix
        host_config_path = os.path.abspath(join_path(dest_dir, host_config_path))

        return host_config_path

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

        if "+openmp" in spec:
            if spec.satisfies("%xl"):
                entries.append(cmake_cache_string("TETON_OpenMP_Fortran_FLAGS_RELEASE", "-qsmp=omp"))
                entries.append(cmake_cache_string("TETON_OpenMP_Fortran_FLAGS_DEBUG", "-qsmp=omp"))
                entries.append(cmake_cache_string("TETON_OpenMP_CXX_FLAGS_RELEASE", "-qsmp=omp"))
                entries.append(cmake_cache_string("TETON_OpenMP_CXX_FLAGS_DEBUG", "-qsmp=omp"))

                if "+openmp" in spec:
                    entries.append(cmake_cache_string("TETON_OpenMP_Offload_Fortran_FLAGS", "-qoffload -qtgtarch=sm_70"))
            else:
                entries.append(cmake_cache_string("TETON_OpenMP_Fortran_FLAGS_RELEASE", "-fopenmp"))
                entries.append(cmake_cache_string("TETON_OpenMP_Fortran_FLAGS_DEBUG", "-fopenmp"))
                entries.append(cmake_cache_string("TETON_OpenMP_CXX_FLAGS_RELEASE", "-fopenmp"))
                entries.append(cmake_cache_string("TETON_OpenMP_CXX_FLAGS_DEBUG", "-fopenmp"))

        if spec.satisfies("%cce"):
            entries.append(cmake_cache_option("STRICT_FPP_MODE", True))
            entries.append(cmake_cache_string("TETON_OpenMP_Offload_LINK_FLAGS", "-ldl"))

        return entries

    def initconfig_hardware_entries(self):
        spec = self.spec
        entries = super(Umt, self).initconfig_hardware_entries()

        #######################
        # Parallelism
        #######################
        if "+cuda" in spec:
            entries.append(cmake_cache_option("ENABLE_CUDA", True))
        if "+openmp" in spec:
            entries.append(cmake_cache_option("ENABLE_OPENMP", True))
        if "+openmp_offload" in spec:
            entries.append(cmake_cache_option("ENABLE_OPENMP_OFFLOAD", True))

        if spec.satisfies("%cce"):
            entries.append(cmake_cache_option("OPENMP_HAS_USE_DEVICE_ADDR", True))
            entries.append(cmake_cache_option("OPENMP_HAS_FORTRAN_INTERFACE", True))

        return entries

    def initconfig_mpi_entries(self):
        entries = super(Umt, self).initconfig_mpi_entries()

        return entries

    def initconfig_package_entries(self):
        spec = self.spec
        entries = []

        #######################
        # Disable features not needed by UMT
        #######################
        entries.append(cmake_cache_option("ENABLE_MINIAPP_BUILD", True))
        entries.append(cmake_cache_option("ENABLE_TESTS", True))
        entries.append(cmake_cache_option("ENABLE_SILO", False))
        entries.append(cmake_cache_option("ENABLE_HDF5", False))

        entries.append(cmake_cache_path("CONDUIT_ROOT", spec["conduit"].prefix))

        if "+mfem" in spec:
            entries.append(cmake_cache_option("ENABLE_MFEM", True))
            entries.append(cmake_cache_path("MFEM_ROOT", spec["mfem"].prefix))
            if "hypre" in spec:
                entries.append(cmake_cache_path("HYPRE_ROOT", spec["hypre"].prefix))
            if "metis" in spec:
                entries.append(cmake_cache_path("METIS_ROOT", spec["metis"].prefix))
            
        if "+caliper" in spec:
            entries.append(cmake_cache_option("ENABLE_CALIPER", True))
            entries.append(cmake_cache_path("CALIPER_ROOT", spec["caliper"].prefix))
            if "adiak" in spec:
                entries.append(cmake_cache_path("ADIAK_ROOT", spec["adiak"].prefix))

        if "+umpire" in spec:
            entries.append(cmake_cache_option("ENABLE_UMPIRE", True))
            entries.append(cmake_cache_path("UMPIRE_ROOT", spec["umpire"].prefix))
            if "camp" in spec:
                entries.append(cmake_cache_option("ENABLE_CAMP", True))
                entries.append(cmake_cache_path("CAMP_ROOT", spec["camp"].prefix))

        return entries
