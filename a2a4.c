#include "a2a4.h"

int main(int argc, char *argv[]) {

  char str[256];

  FILE *dataFile;

  if (argc != 2) {
    printf("\npeak_comp_poisson parameter_file\n");
    printf("Calculates a2 and a4 for angular distribution data.\n\nData format is (in plaintext, one line per value):\nAngle(degrees)   Distribution value   Distribution value error\n\n");
    exit(-1);
  }

  numDataPts=0;
  if ((dataFile = fopen(argv[1], "r")) == NULL) {
    printf("ERROR: Cannot open the data file: %s\n", argv[1]);
    exit(-1);
  }
  while (!(feof(dataFile))) // go until the end of file is reached
  {
    if (fgets(str, 256, dataFile) != NULL) {
        // import angular distribution data
        if (sscanf(str, "%lf %lf %lf", &angle[numDataPts], &val[numDataPts], &err[numDataPts]) == 3) {
          //printf("here!\n");
          cosangle[numDataPts]=cos(angle[numDataPts]*3.14159265/180.0);
          //printf("here!\n");
          numDataPts++;
        }
    }
  }

  printf("Data read in.\n");
  printf("Angle     CosAngle Val      Err\n");
  for(int i=0;i<numDataPts;i++){
    printf("%f %f %f %f\n",angle[i],cosangle[i],val[i],err[i]);
  }
  printf("\n");
  
  find_chisqMin();

  return 0; // great success
}

double a2a4(const double *par) {
  // neyman chisq
  // for more information see Baker and Cousins pg. 438 and Appendix
  double yi = 0.;     // model
  double ni = 0.;     // experiment
  double chisq = 0.; // neyman ratio chisq
  int i = 0;

  for (i = 0; i < numDataPts; i++) {
    // value in ith bin
    ni = val[i];
    // calculate model in the ith bin
    yi = 1.0 + par[0]*0.5*(3.0*cosangle[i]*cosangle[i] - 1.0) + par[1]*0.125*(35.0*cosangle[i]*cosangle[i]*cosangle[i]*cosangle[i] - 30.0*cosangle[i]*cosangle[i] + 3.0);
    yi = yi*par[2];
    //printf("n[%i]: %f, y[%i]: %f, err[%i]: %f\n",i,ni,i,yi,i,err[i]);

    // evaluate chisq given input parameters
    chisq += (ni - yi) * (ni - yi) / (err[i]*err[i]);
    
  }
  //printf("chisq: %f\n",chisq);
  return chisq;
}

void find_chisqMin() {

  printf("Fitting...\n");

  // for more information see minimizer class documentation
  // https://root.cern.ch/root/html/ROOT__Math__Minimizer.html
  char minName[132] = "Minuit";
  char algoName[132] = ""; // default (Migard for Minuit)
  ROOT::Math::Minimizer *min =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    // set tolerance , etc...
    min->SetMaxFunctionCalls(10000000); // for Minuit
    min->SetMaxIterations(100000);
    min->SetTolerance(0.01);
    //min->SetPrintLevel(1); //set to 0 for less info
    min->SetPrintLevel(0); // set to 1 for more info

    // create function wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor lr(&a2a4,3); // neyman chisq */

    // step size and starting variable values
    double step[3] = {0.01,0.01,0.1};
    double variable[3] = {0.0,0.0,10.0};

    min->SetFunction(lr);

    // Set pars for minimization
    min->SetVariable(0, "a2", variable[0], step[0]);
    min->SetVariable(1, "a4", variable[1], step[1]);
    min->SetVariable(2, "scale", variable[2], step[2]);

    // variable limits (optional)
    /* min->SetVariableLimits(0,0.,1E3); */
    /* min->SetVariableLimits(1,0.,1E3); */
    /* min->SetVariableLimits(2,0.,1E3); */

    // do the minimization
    min->Minimize();

    // grab parameters and parameter errors from minimum
    const double *xs = min->X();
    const double *exs = min->Errors();

    // print results
    // degrees of freedom assuming 2 pars
    double ndf = numDataPts - 3;
    
    printf("a2: %f +/- %f\n",xs[0],exs[0]);
    printf("a4: %f +/- %f\n",xs[1],exs[1]);
    printf("scaling: %f +/- %f\n",xs[2],exs[2]);
    printf("chisq: %f\n",min->MinValue());
    printf("chisq/ndf: %f\n",min->MinValue()/ndf);

}