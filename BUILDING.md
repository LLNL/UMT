The Spack package manager has a UMT package and can be used to build UMT and all its required libraries.
For more info on spack see https://github.com/spack/spack

Alternatively, if you want to build UMT and its required libraries by hand, see the script 'build_and_run_umt.sh'.  It is recommended that you are conversant with CMake.

Brief walkthrough on building with Spack.

1. git clone spack from github
``` git clone https://github.com/spack/spack.git ```

2. Setup your environment for spack.  Spack provides shell scripts for bash, csh, 
``` source spack/share/spack/setup-env.sh ```

3. Create an environment for umt
``` spack env create myumt ```

4. Activate your environment
``` spack env activate myumt ```

5. Tell spack to find and add support for your compiler.  The below command will add all compilers that it finds.  Be sure that your compiler is in your path, for spack to find it.  This example assumes that you have gcc version 8.1.0 in your path.
``` spack compiler find```

6. Add your mpi to the spack environment.  This example is using mvapich2 but other mpi versions, such as openmpi, are fine also.
``` spack external find mvapich2

6. Add umt to your environment.
``` spack add umt+mfem %gcc@8.1.0 ```

7. Concretize your environment
``` spack concretize ```

Spack will list all the libraries and dependencies it needs to build.  This list can be reduced by adding more external packages to your spack environment and spack will make use of them instead of building them from scratch.  The perl library is one example that is needed and frequently available on systems already.

```
[me@mymachine:spack]$ spack external find perl
==> The following specs have been detected on this system and added to /home/me/git/spack/var/spack/environments/myumt/spack.yaml
perl@5.16.3
```

Try concretizing again.
``` spack concretize -f ```

8. Build umt
```spack install -j NN ```
where NN is the number of make tasks you want to use.
