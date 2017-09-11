Code example for GEOPHYSICS manuscript GEO-2015-0220  
"CIDRe - a parameter constrained irregular resampling method for scattered point data"

author:       Dipl.-Geoinf. Peter Menzel 
institution:  Christian-Albrechts-Universität zu Kiel, Germany,
              Department for Geosciences

Before using this software, you must read and accept the terms in:
SEG_LEGAL_DISCLAIMER

The disclaimer may also be found at:
http://software.seg.org/disclaimer.txt

This software may be found at:
http://software.seg.org/2016/0002

This example reproduces the results shown in the manuscript for the given synthetic example,
to show how CIDRe is implemented and has to be parameterized. Additionally, it demonstrates
in which way the other algorithms are used for comparison. Data loading, parameterizing and
call of all algorithms used and the comparison of the results are performed in MainScript.m

This code consists of 3 code files:

MainScript.m:               Main file to be executed in Matlab.
			    Here the data are loaded and CIDRe is parameterized.
CIDRe2D.m:                  Methods for CIDRe algorithm.
Other_Resampling_Methods.m: Other methods used for comparison.

and 1 data file:

syntheticExample.csv:       Synthetic data in comma separated format [x,y,value].

All 4 files have to be located in the same folder.
The code is implemented in Matlab; it has to be executed within the Matlab software.
To execute the main script, change the active Matlab directory to the containing
folder and type:

>>MainScript

into the Matlab command window or open the file in the Matlab-Editor and run it (F5). 
After execution, the script produces the figures 5 a,b,c,d, described in the manuscript 
(Section 'VALIDATION BY SYNTHETIC DATA'). The parameters used to create the results shown 
in for Figure 4 are given in MainScript.m as comments.

Note: The code is not fully optimized and not parallelized yet.
For the complete script an over all calculation time of <2000s will be needed
on common hardware (circa 2015).
