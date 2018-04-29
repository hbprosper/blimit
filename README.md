# blimit


## Description 

A simple Bayes limit calculator for the single count model.
The likelihood function is assumed to be

      Pr(n|mean) = exp(-mean) mean^n/n!

where

    mean = a * x + b and "x" is the cross-section.

__Note__: if *a* is set to the expected signal, then *x* will be the
signal strength

	r = signal / predicted signal.

The prior for *a* is assumed to be a gamma density derived from
the error and estimate of *a*. Likewise for *b*. We assume that
effective count in the gamma density is given by

	Neff = (estimate/error)^2.

This follows from the assumption that the estimate of *a* is derived by
scaling down the count *Neff* by the factor *k = estimate / Neff*.

*blimit* can be used for problems in
which Gaussian priors are inappropriate. This is typically the case when the
Gaussian would be severely truncated at one bound of the parameter space.
*blimit* uses gamma priors for the effective integrated luminosity (defined
as the acceptance * integrated luminosity) and for the background, namely, a
prior that goes to zero at zero effective luminosity and zero background. 
Compute cross-section limit given observed count n, estimates for the
acceptance * times luminosity, *a*, and background, *b*, with known errors. 
The prior for the cross-section is taken to be flat up to some upper bound.

 
## Examples

Build the program blimit using

	make

and add to PATH using

	source setup.sh

Usage:

        blimit  -n <count>
	        -l <effective lum[,error]>
                -b [background (0)[,error]]
                -c [CL (0.95)]

If the background is not given, it is assumed to be zero (actually 1.e-4). 

### Example 1:

	blimit -n0 -l1 -c0.9 

	observed count:        0
	eff.lumi/sigmal:       1.000 +/- 0.000
	background:            0.000 +/- 0.000
	CL:                    0.900
	Searching for limit in range [0, 8]
	upper limit:        2.2996

as one expects for zero background, unit effective integrated
luminosity (or signal strength), and
zero relative error and a flat prior for the cross-section.

### Example 2:

	blimit -n1 -l1 -c0.90

	observed count:        1
	eff.lumi/sigmal:       1.000 +/- 0.000
	background:            0.000 +/- 0.000
	CL:                    0.900
	searching for limit in range [0, 12.3137]
	upper limit:        3.8889


### Example 3:

	blimit -n1 -l1,0.2 -c0.90

	observed count:        1
	eff.lumi/sigmal:       1.000 +/- 0.200
	background:            0.000 +/- 0.000
	CL:                    0.900
	searching for limit in range [0, 12.3137]
	upper limit:        4.1162
