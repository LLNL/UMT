Example of using this spack environment to build on a Cray system with MI300 GPUs with the Caliper timer library and Umpire memory pool library enabled.

1. wget https://github.com/spack/spack/releases/download/v0.23.1/spack-0.23.1.tar.gz
2. tar xvzf spack-0.23.1.tar.gz 
3. cp -r umt ../../spack-0.23.1/var/spack/repos/builtin/packages
4. source spack-0.23.1/share/spack/setup-env.sh
5. spack -e umt/spack/environment/cray install -j40
