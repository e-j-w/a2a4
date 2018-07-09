#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// ROOT libraries
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#define NPTS        10000

//forward declarations
void   find_chisqMin(); // ROOT minimization class

//global variables
double angle[NPTS],cosangle[NPTS],val[NPTS],err[NPTS];
int numDataPts;

// ROOT stuff
double a2a4(const double *par); // function to get chisq
