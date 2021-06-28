# Getting started

_This will build all dependencies required to build OpenFAST, OpenFOAM, and SOWFA._

1. Clone https://github.com/spack/spack

2. Add to `.bash_profile` (or something similar):
```
export SPACK_ROOT=${HOME}/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
```

3. Clone https://github.com/jrood-nrel/spack-configs and run `scripts/setup-spack.sh` to grab `packages.yaml` and `config.yaml` (Note: this depends on `${SPACK_ROOT}` in your environment)

4. Adapt the following commands:
```bash
# install dependencies for OpenFOAM-6
nice spack install --only dependencies openfoam-org@6 +int64 +metis %intel@18.0.4
# install parmetis
nice spack install parmetis@4.0.3 +int64 +shared %intel@18.0.4
# install CGAL (includes dependencies: boost, gmp, mpfr)
nice spack install cgal@4.13 %intel@18.0.4
# install OpenFAST dependencies
nice spack install hdf5@1.10.6 +cxx +fortran +hl -mpi +pic +shared +szip %intel@18.0.4
nice spack install libxml2@2.9.10 +python %intel@18.0.4
nice spack install yaml-cpp@0.6.3 %intel@18.0.4
nice spack install zlib@1.2.11 %intel@18.0.4
```
You should make sure there are no errors from the installs.
