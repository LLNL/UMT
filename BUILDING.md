The Spack package manager has a UMT package and can be used to build UMT and all its required libraries.
For more info on spack see https://github.com/spack/spack

Alternatively, if you want to build UMT and its required libraries by hand, see the script 'build_umt.sh'.  It is recommended that you are conversant with CMake.


Brief walkthrough on building with Spack.

1. git clone spack from github
``` https://github.com/spack/spack.git ```

2. Setup your environment for spack
``` source spack/var/spack/setup-env.sh ```

3. Create an environment for umt
``` spack env create myumt ```

4. Activate your environment
``` spack env activate myumt ```

5. Depending on your compiler, spack may have it already found.  If not, you may need to add it manually.  See the spack documentation for more details on this.
``` spack compiler find ```

6. Add umt to your environment.  This examples is using gcc 8.1.0
``` spack add umt+mfem %gcc@8.1.0 ```

7. Concretize your environment
``` spack concretize ```

Spack will list all the libraries and dependencies it needs to build.  This may be substantial if Spack is not aware of common libraries your platform already has available.  To tell Spack to try and find these present libraries, use 'spack external find'.  For example, these are common libraries that may be present on your platform already:

```
[me@mymachine:spack]$ spack external find ncurses
==> The following specs have been detected on this system and added to /home/me/git/spack/var/spack/environments/myumt/spack.yaml
ncurses@5.9.20130511
[me@mymachine:spack]$ spack external find openssl
==> The following specs have been detected on this system and added to /home/me/git/spack/var/spack/environments/myumt/spack.yaml
openssl@1.0.2k-fips
[me@mymachine:spack]$ spack external find mvapich2
==> The following specs have been detected on this system and added to /home/me/git/spack/var/spack/environments/myumt/spack.yaml
mvapich2@2.3.1
[me@mymachine:spack]$ spack external find perl
==> The following specs have been detected on this system and added to /home/me/git/spack/var/spack/environments/myumt/spack.yaml
perl@5.16.3
```

Try concretizing again.
``` spack concretize -f ```

8. Build umt
```spack install -j NN ```
where NN is the number of make tasks you want to use.
