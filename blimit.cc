////////////////////////////////////////////////////////////////////////////
// File: blimit.cc
//
// Description: 
//
//  Compute cross-section limit given observed count n and estimates for the
//  acceptance * times luminosity, a, and background, b. The prior for the 
//  cross-section is taken to be flat.
//
// Notes:
// (A)
//    The likelihood function is assumed to be
//
//      Pr(n|mean) = exp(-mean) mean^n/n!
//
//    where
//
//    mean = a * x + b and "x" is the cross-section.
//
// (B)
//    The prior for a is assumed to be a gamma density derived from
//    the error and estimate of a. Likewise for b. We assume that the 
//    effective count is given by
//
//        A = (estimate/error)^2
//
//    and the scale factor that scales down the effective count is given by
//
//        alpha = estimate/A = error^2 / estimate
//
// Created  6-Jun-2000  Harrison B. Prosper
// Updated  6-Dec-2000  HBP, Make separate "eps" value for integration and
//                      root finding and use double precision where possible
// Updated  6-Jun-2001  HBP, Make it possible to use non-integral numbers of
//                      "observed" events.
// Updated  6-Jun-2005  HBP this is basically a simplified version of climit
//                      but which is valid for all error sizes.
//         15-Jun-2005  HBP isolate main program. 
//         19-Oct-2005  HBP add some checks.
//         28-Nov-2005  HBP compute beta functions another way
//         02-Mar-2011  HBP use MathCore classes. Number of observed events
//                      must be integer. Go back to PoissonGamma calculation
////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "blimit.h"

#include "TMath.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
#include "Math/RootFinder.h"

using namespace std;
//--------------------------------------------------------------------------
// Hide in anonymous namespace
//--------------------------------------------------------------------------
namespace 
{
  const  int   MAXCOUNT  = 10000;
  const  float DEFCL     = 0.95;     // Default confidence level
  const  double ABSTOL   = 1.e-9;
  const  double RELTOL   = 1.e-6;
  const  int    SIZE     = 100;
  const  int    RULE     = 3;

  struct Post
  {
    int    D;
    double CL;                 // Confidence Level
    double XMIN;
    double XMAX;
    double NORM;

    int N;
    double eA;
    double Aeff;

    double eB;
    double Beff;

    vector<double> p;
    vector<double> y;
    vector<double> A;
    vector<vector<long double> > c;

    ROOT::Math::WrappedMemFunction<Post, double (Post::*)(double)> fn;
    ROOT::Math::IntegratorOneDim ifn;

    Post(int count, 
         double efflum, double efflumError,
         double bkg, double bkgError,
         double cl, 
         double xmin, double xmax)
      : D(count),
        CL(cl),
        XMIN(xmin),
        XMAX(xmax),
        NORM(1),
        fn(*this, &Post::density),
        ifn(fn, 
            ROOT::Math::IntegrationOneDim::kADAPTIVE,
            ABSTOL,
            RELTOL,
            SIZE,
            RULE)
    {
      // Compute scale factors and effective counts
      A.clear();
      p.clear();
      y.clear();
      c.clear();

      if ( efflum <= 0 )
        {
          efflum = 1;
          efflumError = 1.e-4 * efflum;
        }
      if ( efflumError <= 0 ) efflumError = 1.e-4 * efflum;
      Aeff = efflum / efflumError;
      Aeff*= Aeff;
      eA   = efflum / Aeff;
      A.push_back(Aeff);
      
      N = 1;
      p.push_back(0);
      y.push_back(0);
      c.push_back(vector<long double>(D+1));

      if ( bkg <= 0 )
        {
          bkg = 1.e-4;
          bkgError = 1.e-4 * bkg;
        }
      if ( bkgError <= 0 ) bkgError = 1.e-4 * bkg;
      Beff = bkg / bkgError;
      Beff*= Beff;
      eB   = bkg / Beff;
      A.push_back(Beff);
      
      N++;
      p.push_back(eB);
      y.push_back(0);
      c.push_back(vector<long double>(D+1));

//       printf("efflum = %10.3e\tefflumError = %10.3e\n"
//              "bkg    = %10.3e\tbkgError    = %10.3e\n", 
//              efflum, efflumError, bkg, bkgError);

//       printf("Aeff = %10.3e\teA = %10.3e\n"
//              "Beff = %10.3e\teB = %10.3e\n", Aeff, eA, Beff, eB);

      // compute normalization
      
      NORM = 1;
      NORM = ifn.Integral(XMIN, XMAX);

//       printf("Norm = %10.3e\n", NORM);
    }

    ~Post() {}

    //----------------------------------------------------------------------
    // Posterior Density Function
    //----------------------------------------------------------------------
    
    double density(double xsec)
    {
      p[0] = xsec * eA;
      p[1] = eB;

      for (int j = 0; j < N; ++j) y[j] = p[j] / (p[j] + 1);

      // compute coefficients
      
      for (int j = 0; j < N; ++j) c[j][0] = pow(1 + p[j], -(A[j] + 1));
      
      for (int k = 1; k < D+1; ++k)
        for (int j = 0; j < N; ++j)
          c[j][k] = c[j][k-1] * y[j] * (A[j] + k) / k;

      // compute prob by summing D+1 terms
      double prob = 0.0;
      for (int j = 0; j <= D; ++j) prob += c[0][j] * c[1][D-j];

      return prob;
    }

    double cdf(double x) { return ifn.Integral(XMIN, x) / NORM; } 

    double limit()
    {
      // function whose root is to be found
      ROOT::Math::WrappedMemFunction<Post,
        double (Post::*)(double)> fn(*this, &Post::f);
      ROOT::Math::RootFinder rootfinder;
  
        cout << "searching for limit in range [" 
           << XMIN << ", " 
           << XMAX << "]" << endl;

      rootfinder.SetFunction(fn, XMIN, XMAX);
      int status = rootfinder.Solve();
      if ( status != 1 )
      {
        cout << "*** Post *** RootFinder failed"
             << endl;
        return -1;
      }
      return rootfinder.Root();
    }

    //----------------------------------------------------------------------
    // Function whose root is to be found f(x) = CDF(x) - CL = 0
    //----------------------------------------------------------------------
    double f(double x)
    {
      return cdf(x) - CL;
    }
  };

  char replace(char s) 
  { 
    if ( s == ',' ) 
      return ' ';
    else
      return s;
  }

  void decodet(string line, float& estimate, float& error)
  {
    //cout << line << endl;
    estimate = 0;
    error = 0;
    transform(line.begin(), line.end(), line.begin(), replace);
    istringstream stream(line.c_str());
    stream >> estimate >> error;
  }
};

//--------------------------------------------------------------------------
// Exported interface
//--------------------------------------------------------------------------
float blimit(int   count,
             float efflum,
             float efflumError,
             float bkg,
             float bkgError,
             float cl)
{

  // Check for sensible count

  if ( count < 0 )
    {
      cout << "Observed count should be positive!" << endl;
      return -1;
    }
  if ( count > MAXCOUNT )
    {
      cout << "Observed count (" << count << ") too large!" << endl;
      return -1;
    }


  // Check effective luminosity

  if ( efflum < 0 )
    {
      cout << "Effective integrated luminosity should be positive!" << endl;
      return -1;
    }
  
  // Check background

  if ( bkg < 0 )
    {
      cout << "Background should be positive!" << endl;
      return -1;
    }

  if ( cl < 0 || cl > 1 )
    {
      cout << "CL must be between 0 and 1!" << endl;
      return -1;
    }

  // Compute upper limit

  float sigma = sqrt((double)(count + bkgError*bkgError + 1));
  float xmin  = 0;
  float xmax  = (fabs(count - bkg) + 8*sigma)/efflum;

  Post post(count,
            efflum, efflumError,
            bkg, bkgError,
            cl,
            xmin, xmax);
  return post.limit();
}

#ifdef __MAIN__
//--------------------------------------------------------------------------
// MAIN ROUTINE
//--------------------------------------------------------------------------
int main(int argc, char** argv)
{
  //---------------------------------------------
  // Define usage string
  //---------------------------------------------
  string usage("Usage: \n  ./blimit");
  usage += "\t-n <count>\n";
  usage += "\t\t-l <effective lum[,error]> (1)\n";
  usage += "\t\t-b <background (0)[,error]> (0)\n";
  usage += "\t\t-c <CL> (0.95)\n";

  if ( argc < 1 )
    {
      cout << usage << endl;
      return -1;
    }

  //---------------------------------------------
  // Process inputs
  //---------------------------------------------
  int option;

  map<string, string> opt;
  double estimate, error;
  while ( (option = getopt(argc, argv, "l:b:n:c:h")) != EOF )
    {
      switch (option)
        {
        case 'l': opt["lum"]  = string(optarg); break;
        case 'b': opt["bkg"]  = string(optarg); break;
        case 'n': opt["count"]= string(optarg); break;
        case 'c': opt["CL"]   = string(optarg); break;
        case 'h': cout << usage;            return 0;
        }
    }

  //---------------------------------------------
  // Check arguments and set defaults
  //---------------------------------------------

  int   count;
  float efflum=1;
  float efflumError=0;
  float bkg=0;
  float bkgError=0;
  float cl=0.95;

  // Check observed count

  if ( opt.find("count") == opt.end() )
    {
      cout << "*** Need observed count"
	   << endl << usage << endl;
      return -1;
    }
  else
    count = atoi(opt["count"].c_str());
  if ( count < 0 )
    {
      cout << "*** Observed count must >= 0" << endl;
      return -1;
    }

  // Check effective luminosity

  if ( opt.find("lum") != opt.end() )
    decodet(opt["lum"], efflum, efflumError);
  if ( efflum <= 0.0 )
    {
      cout << "*** Effective integrated luminosity must > 0" << endl;
      return -1;
    }

  // Check background

  if ( opt.find("bkg") != opt.end() )
    decodet(opt["bkg"], bkg, bkgError);
  if ( bkg < 0.0 )
    {
      cout << "Background must >= 0" << endl;
      return -1;
    }


  if ( opt.find("CL") != opt.end() )
    cl = atof(opt["CL"].c_str());
  else
    cl = DEFCL;

  char str[800];
  sprintf(str,
          "observed count:  %7d\n"
          "eff.lumi/sigmal:  %10.3f +/- %5.3f\n"
          "background:       %10.3f +/- %5.3f\n"
          "CL:               %10.3f",
          count,
          efflum, efflumError,
          bkg, bkgError,
          cl);
  
  cout << endl << str << endl;

  // Compute upper limit

  double xsecup= blimit(count,
                        efflum, efflumError,
                        bkg, bkgError,
                        cl);
  // Result
  
  sprintf(str, "upper limit:    %10.4f\n", xsecup);
  cout << str << endl;

  return 0;  
}
#endif
