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
#define BIG_NUMBER  1E30

//forward declarations
void   find_chisqMin(int); // ROOT minimization class

//global variables
double angle[NPTS],cosangle[NPTS],val[NPTS],err[NPTS];
int numDataPts;
double fixedA2,fixedA4,fixedA6;

// ROOT stuff
double a2a4(const double *par); // function to get chisq

enum fit_mode_enum{FITMODE_A2, FITMODE_A4, FITMODE_A6, FITMODE_A2A4, FITMODE_A2A4A6, FITMODE_A2NOSCALING, FITMODE_ENUM_LENGTH};
