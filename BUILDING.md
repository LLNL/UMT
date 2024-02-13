UMT requires CMake to build.  It is recommended that your are familiar with CMake, in order to enable/disable the supported build options for UMT.  An example shell script 'build_and_run_umt.sh' is included in the UMT that will compile UMT, and its required libraries, on a typical Linux distribution.  The Linux distribution must have an MPI installation available.

A list of common CMake options for UMT, and their default values, are:

Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel
**CMAKE_BUILD_TYPE:STRING=**

Link and enable the use of the optional CALIPER performance measurement library.
**ENABLE_CALIPER:BOOL=OFF**

Enable compiling CUDA code
**ENABLE_CUDA:BOOL=OFF**

Enables CMake's Find MPI support (Turn off when compiling with the mpi wrapper directly)
**ENABLE_FIND_MPI:BOOL=ON**

Enable OpenMP pragmas
**ENABLE_OPENMP:BOOL=OFF**

Enable OpenMP target offload pragmas
**ENABLE_OPENMP_OFFLOAD:BOOL=OFF**

Link and enable the use of the optional UMPIRE memory library.  Umpire provides memory pools on both CPU and GPU devices, and CPU memory pools in page-locked memory.
**ENABLE_UMPIRE:BOOL=OFF**

Compiler OpenMP runtime implementation includes Fortran interface
**OPENMP_HAS_FORTRAN_INTERFACE:BOOL=OFF**

Compiler OpenMP implementation includes use_device_addr pragma.  This is usually true for OpenMP 5.x implementations, but not for OpenMP 4.x.
**OPENMP_HAS_USE_DEVICE_ADDR:BOOL=OFF**

Use preprocessor macros that strictly conform to FPP.
**STRICT_FPP_MODE:BOOL=OFF**





UMT also provides a Spack package and can be optionally built using that package manager.  For more information on Spack see https://github.com/spack/spack.

For required dependencies on UMT, see DEPENDENCIES.md.

BUILDING WITH SHELL SCRIPT:
1. Open the bash script 'build_and_run_umt.sh' and refer to the instruction there-in.

If you are testing code changes to UMT, you can make needed changes and recompile the code by running 'make install' in the build directory.

BUILDING WITH SPACK:
It is recommended that developers are familiar with Spack before trying to build UMT with this package manager.  Spack provides a tutorial at https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html.  A brief walkthrough is provided here for convenience.

1. git clone spack from github
``` git clone https://github.com/spack/spack.git ```

2. Setup your environment for spack.  Spack provides shell scripts for bash, csh, 
``` source spack/share/spack/setup-env.sh ```

3. Create an environment for umt
``` spack env create myumt ```

4. Activate your environment
``` spack env activate myumt ```

5. Tell spack to find and add support for your compiler.  The below command will add all compilers that it finds.  If Spack does not find your compiler, try adding it to your path.
``` spack compiler find```

6. Add your mpi to the spack environment.  This example is using mvapich2.
``` spack external find mvapich2 ```

7. Add cmake to your spack environment.  Note: UMT requires at least CMake 3.14, so be sure Spack finds 3.14 or newer.
``` spack external find cmake ```

8. Add umt to your environment.
``` spack add umt+mfem %gcc@8.1.0 ^mvapich2```

9. Concretize your environment
``` spack concretize ```

Spack will list all the libraries and dependencies it needs to build.  This list can be reduced by adding more external packages to your spack environment and spack will make use of them instead of building them from scratch.  The perl library is one example that is needed by the Hypre library and is frequently available on systems already.  An example of adding perl, then concretizing again, is below.

```
[me@mymachine:spack]$ spack external find perl
==> The following specs have been detected on this system and added to /home/me/git/spack/var/spack/environments/myumt/spack.yaml
perl@5.16.3
[me@mymachine:spack]$ spack concretize -f
```

8. Build umt
```spack install -j NN ```
where NN is the number of make tasks you want to use.

UMT will now be installed in your spack installation directory under "spack/opt/spack/<platform>/<compiler>"

If you are testing code changes to UMT, you can make needed changes and recompile the code using Spack by using its 'developer' mode.  See the Spack tutorial on this feature at https://spack-tutorial.readthedocs.io/en/latest/tutorial_developer_workflows.html.
