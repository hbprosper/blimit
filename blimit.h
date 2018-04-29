#ifndef BLIMIT_H
#define BLIMIT_H
////////////////////////////////////////////////////////////////////////////
// File: blimit.h
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
//    scale factor is given by
//
//        e = sigma^2 / estimate
//
// Created  6-Jun-2000  Harrison B. Prosper
// Updated  6-Dec-2000  HBP, Make separate "eps" value for integration and
//                      root finding and use double precision where possible
// Updated  6-Jun-2001  HBP, Make it possible to use non-integral numbers of
//                      "observed" events.
// Updated  6-Jun-2005  HBP this is basically a simplified version of climit
//                      but which is valid for all error sizes.
//         15-Jun-2005  HBP isolate main program. 
////////////////////////////////////////////////////////////////////////////

float blimit(int   count,
             float efflum,
             float efflumError,
             float bkg,
             float bkgError,
             float cl=0.95);

#endif
