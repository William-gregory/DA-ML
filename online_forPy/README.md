Online implementation of an ML-based correction to category sea ice concentrations in SIS2; implemented via Forpy. This repository contains the Python script `NNetwork` which is called in Fortran. The repositories which contain the Fortran code to execute this are in [SIS2](https://github.com/William-gregory/SIS2/tree/forpy_dev/src).

Specifically, see source code:

- [SIS_slow_thermo.F90](https://github.com/William-gregory/SIS2/tree/forpy_dev/src/SIS_slow_thermo.F90)
- [SIS_G23_CNN.F90](https://github.com/William-gregory/SIS2/tree/forpy_dev/src/SIS_G23_CNN.F90)
- [Forpy_interface.F90](https://github.com/William-gregory/SIS2/tree/forpy_dev/src/Forpy_interface.F90)
- [forpy_mod.F90](https://github.com/William-gregory/SIS2/tree/forpy_dev/src/forpy_mod.F90)

For compiling the MOM6/SIS2 model code in order to run Forpy, see the steps outlined in [MOM6-examples](https://github.com/William-gregory/MOM6-examples/blob/forpy_dev/COMPILE_MOM6SIS2.sh).
