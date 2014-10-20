# Fast Asimov Utils

*Kyle Cranmer & Glen Cowan*

*License LGPL v2.1 (for ROOT compatibility, happy to make it BSD for other purposes)*


This code quickly calculates the expected 95% CLs upper limit on the number of
signal events, `s`, given an expected background, `bExp`, and 
uncertainty on the background estimate, `deltaB`.
The background uncertainty is absolute (not relative) and is uncertainty
on the mean background (so you don't include sqrt(bExp) Poisson fluctuations in this number).

Example Usage: 

  * you expect 50 +/- 7 background events
  * `ExpectedLimit(50,7)`
  * returns s_95 = 19.7 events

Similar code for expected discovery significance is also included.

Example Usage:

   * you expect 50 signal events, 100 +/- 7 background events
   * `ExpectedSignificance(50,100,7)`
   * returns 3.72 sigma


The derivations of these formulae are based on a statistical model:

	Pois(n | s+b ) * Pois(m | tau * b)

The tau quantity is calculated from `tau=bExp/deltaB/deltaB`.
See arXiv:physics/0702156 for motivation of this model.
The maximum likelihood estimate and conditional maximum likelihood estimate
were solved analytically and coded here.
The log-likelihood ratio and the profile log likelihood ratio follow immediately.
The b-only and s+b p-values needed for CLs are calculated using the 
asymptotic distributions in 
>	Cowan, Cranmer, Gross, Vitells,	Eur. Phys. J. C 71 (2011) 1554.  
[arXiv:1007.1727](http://arxiv.org/abs/1007.1727)

The ExpectedSignificance formulae was derived by Cowan and the numerical solution 
for the Expected upper limit was written by Cranmer as correlaries to that paper.

Note, the ExpectedSignificance equation is the same as Eq.(17) of
> Tipei Li and Yuqian Ma, Astrophysical Journal 272 (1983) 317–324.

and Eq.(25) of 
>  Robert D. Cousins, James T. Linnemann and Jordan Tucker, NIM A 595 (2008) 480– 501; [arXiv:physics/0702156](http://arxiv.org/abs/physics/0702156).

after making the replacements `n=bExp` and `m=bExp*tau`. 

## Validation

This code could use more validation, but there are a series of tests that run that check against results from the mathematica notebook and the equivalent HistFactory example.
The HistFactory validation is in the directory `validation` and includes:
   * A HistFactory model in xml that defines the statistical model above
   * The top-level script `makeHists.C` that writes histograms to file `data/histograms.root` based on the values of `bExp` and `deltaB`, then runs HistFactory's `hist2workspace` + `RooStats/StandardHypoTestInvDemo.C` with the `AsymptoticCalculator` using the 1-sided profile likelihood ratio test statistic and CLs=True.


## To Do

This is a work in progress!
   * C++ to be updated to use Brent q root finding instead of simple scan
   * Could remove ROOT dependency in C++ verison entirely
   * add +/- 1,2 sigma bands for the expected upper limit
   * could do same for observed upper limit.

Extension:
   * Include uncertainty on signal efficiency for upper limit. Those equations are big! (See ExpectedLimitWithSignalUncert.nb)
