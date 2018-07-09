# A2A4

Angular Distribution A2/A4 fitter

Maintainer: Jonathan Williams

## Description

Fits A2, A4 angular distribution coefficients to experimental data using Minuit (ROOT).  See the file `sample_data.txt` for example data (angles specified in degrees).

### Useful References: 

* Phys. Rev. C 45, 2161 (1992)
* Phys. Rev. C 89, 024324 (2014)

## Installation

Requires ROOT (tested with v5.x) to be installed with environment variables set up properly.  Environment variable setup can be done by adding to your `.bashrc` (and then reloading the terminal):

```
#ROOT configuration in .bashrc
export ROOTSYS=/path/to/root
export ROOTINC=$ROOTSYS/include
export ROOTLIB=$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTLIB
```

Once this is done, use `make` to compile.  Tested using g++ and GNU make on Ubuntu 16.04.
