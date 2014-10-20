# Fast Asimov Utils

*Kyle Cranmer & Glen Cowan *

*License LGPL v2.1 (for ROOT compatibility, happy to make it BSD for other purposes) *

This code quickly calculates the expected 95% CLs upper limit on the number of
signal events, `s`, given an expected background, `bExp`, and 
uncertainty on the background estimate, `deltaB`.
The background uncertainty is absolute (not relative) and is uncertainty
on the mean background (so you don't include Poisson fluctuatiosn in this number).
Example Usage: 

  * you expect 50 +/- 3 background events
  * `ExpectedLimit(50,3)`
  * returns s_95 = 13.4 events

Similar code for expected discovery significance is also included
Example Usage:

   * you expect 50 signal events, 100 +/- 7 background events
   * `ExpectedSignificance(50,100,7)`
   * returns 3.72 sigma


The derivations of these formulae are based on a statistical model:

	Pois(n | s+b ) * Pois(m | tau * b)

The tau quantity is calculated from `tau=bExp/deltaB/deltaB`.
The maximum likelihood estimate and conditional maximum likelihood estimate
were solved analytically and coded here.
The log-likelihood ratio and the profile log likelihood ratio follow immediately.
The b-only and s+b p-values needed for CLs are calculated using the 
asymptotic distributions in 
>	Cowan, Cranmer, Gross, Vitells,	Eur. Phys. J. C 71 (2011) 1554.  
[arXiv:1007.1727](http://arxiv.org/abs/1007.1727)

The ExpectedSignificance formulae was derived by Cowan and the numerical solution 
for the Expected upper limit was written by Cranmer as correlaries to that paper.