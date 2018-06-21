# Compiling and Running umt2016

## 0. Choose compilers and cuda version. I have been successful with PGI 18.3 and cuda 9.1.76
```
module load cuda/9.1.76
module load pgi/18.3
```
Make sure make.defs paths agree with the version of pgi and cuda you are using.
Also make sure mpicc and mpif90 use the pgi compiler.

## 1 (optional). If using modules does not automatically identify mpif90->pgf90, then set them:

```
source MPIdefs.sh
```

## 2. Then make in this top directory (make ibm timers first). If changes need to be made to compiler flags, they usually go in make.defs

```
make veryclean
make clean
cd ibmtimers
make
cd ..
make
```

## 3. Link by going to Teton directory and making SuOlsonTest target:

```
cd Teton
make SuOlsonTest
```

## 4. Run a test problem: 

### 4.1 See june-runs for acceptence benchmark run scripts and output suitable for SUMMIT from June runs on PEAK.