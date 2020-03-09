#include "a2a4.h"

int main(int argc, char *argv[]) {

  char str[256];
  int mode; //0=a2a4,1=a2

  FILE *dataFile;

  if ((argc != 2)&&(argc != 3)) {
    printf("\na2a4 parameter_file mode\n");
    printf("Calculates a2 and a4 for angular distribution data.\n\nData format is (in plaintext, one line per value):\nAngle(degrees)   Distribution value   Distribution value error\n\n");
    printf("Valid options for 'mode' are:\na2a4 - Fit both A2 and A4.\na2 - Fit a2 only.\na2noscaling - Fit a2 only (no scaling).\nIf no valid mode is given, a2a4 is assumed.\n\n");
    exit(-1);
  }

  if(argc == 3){
    if(strcmp(argv[2],"a2a4")==0){
      mode=0;
    }else if(strcmp(argv[2],"a2")==0){
      mode=1;
    }else if(strcmp(argv[2],"a2noscaling")==0){
      mode=2;
    }else{
      mode=0;//a2a4 is default
    }
  }else{
    mode=0;//a2a4 is default
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
  
  find_chisqMin(mode);

  return 0; // great success
}

double a2a4(const double *par) {
  
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

double a2(const double *par) {
  
  double yi = 0.;     // model
  double ni = 0.;     // experiment
  double chisq = 0.; // neyman ratio chisq
  int i = 0;

  for (i = 0; i < numDataPts; i++) {
    // value in ith bin
    ni = val[i];
    // calculate model in the ith bin
    yi = 1.0 + par[0]*0.5*(3.0*cosangle[i]*cosangle[i] - 1.0);
    yi = yi*par[1];
    //printf("n[%i]: %f, y[%i]: %f, err[%i]: %f\n",i,ni,i,yi,i,err[i]);

    // evaluate chisq given input parameters
    chisq += (ni - yi) * (ni - yi) / (err[i]*err[i]);
    
  }
  //printf("chisq: %f\n",chisq);
  return chisq;
}

double a2_no_scaling(const double *par) {
  
  double yi = 0.;     // model
  double ni = 0.;     // experiment
  double chisq = 0.; // neyman ratio chisq
  int i = 0;

  for (i = 0; i < numDataPts; i++) {
    // value in ith bin
    ni = val[i];
    // calculate model in the ith bin
    yi = 1.0 + par[0]*0.5*(3.0*cosangle[i]*cosangle[i] - 1.0);
    //printf("n[%i]: %f, y[%i]: %f, err[%i]: %f\n",i,ni,i,yi,i,err[i]);

    // evaluate chisq given input parameters
    chisq += (ni - yi) * (ni - yi) / (err[i]*err[i]);
    
  }
  //printf("chisq: %f\n",chisq);
  return chisq;
}

void find_chisqMin(int mode) {

  printf("Fitting...\n");

  // for more information see minimizer class documentation
  // https://root.cern.ch/root/html/ROOT__Math__Minimizer.html
  char minName[132] = "Minuit";
  char algoName[132] = ""; // default (Migard for Minuit)
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

  // set tolerance , etc...
  min->SetMaxFunctionCalls(10000000); // for Minuit
  min->SetMaxIterations(100000);
  min->SetTolerance(0.01);
  //min->SetPrintLevel(1); //set to 0 for less info
  min->SetPrintLevel(0); // set to 1 for more info

  // create function wrapper for minmizer
  // a IMultiGenFunction type
  if(mode==0){

    ROOT::Math::Functor lr(&a2a4,3); //a2a4 function
    min->SetFunction(lr);

    // step size and starting variable values
    double step[3] = {0.01,0.01,0.1};
    double variable[3] = {0.0,0.0,1.0};

    // Set pars for minimization
    min->SetVariable(0, "a2", variable[0], step[0]);
    min->SetVariable(1, "a4", variable[1], step[1]);
    min->SetVariable(2, "scale", variable[2], step[2]);

    // do the minimization
    min->Minimize();

    // grab parameters and parameter errors from minimum
    const double *xs = min->X();
    const double *exs = min->Errors();

    // print results
    double ndf = numDataPts - 3; //degrees of freedom assuming 3 pars
    printf("Fit function: f(x) = 1 + a2*p2(cos (theta)) + a4*p4(cos (theta))\n");
    printf("a2: %f +/- %f\n",xs[0],exs[0]);
    printf("a4: %f +/- %f\n",xs[1],exs[1]);
    printf("scaling: %f +/- %f\n",xs[2],exs[2]);
    printf("chisq: %f\n",min->MinValue());
    printf("chisq/ndf: %f\n",min->MinValue()/ndf);

  }else if(mode==1){

    ROOT::Math::Functor lr(&a2,2); //a2 function
    min->SetFunction(lr);

    // step size and starting variable values
    double step[2] = {0.01,0.01};
    double variable[2] = {0.0,0.0};

    // Set pars for minimization
    min->SetVariable(0, "a2", variable[0], step[0]);
    min->SetVariable(1, "scale", variable[1], step[1]);

    // do the minimization
    min->Minimize();

    // grab parameters and parameter errors from minimum
    const double *xs = min->X();
    const double *exs = min->Errors();

    // print results
    double ndf = numDataPts - 2; //degrees of freedom assuming 2 pars
    printf("Fit function: f(x) = 1 + a2*p2(cos (theta))\n");
    printf("a2: %f +/- %f\n",xs[0],exs[0]);
    printf("scaling: %f +/- %f\n",xs[1],exs[1]);
    printf("chisq: %f\n",min->MinValue());
    printf("chisq/ndf: %f\n",min->MinValue()/ndf);

  }else if(mode==2){

    ROOT::Math::Functor lr(&a2_no_scaling,1); //a2 function (no scaling)
    min->SetFunction(lr);

    // step size and starting variable values
    double step[1] = {0.01};
    double variable[1] = {0.0};

    // Set pars for minimization
    min->SetVariable(0, "a2", variable[0], step[0]);

    // do the minimization
    min->Minimize();

    // grab parameters and parameter errors from minimum
    const double *xs = min->X();
    const double *exs = min->Errors();

    // print results
    double ndf = numDataPts - 2; //degrees of freedom assuming 2 pars
    printf("Fit function: f(x) = 1 + a2*p2(cos (theta))\n");
    printf("a2: %f +/- %f\n",xs[0],exs[0]);
    printf("chisq: %f\n",min->MinValue());
    printf("chisq/ndf: %f\n",min->MinValue()/ndf);

  }
  
  

}