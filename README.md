# PREAMBLE

This script performs time-domain RFI clipping on 8-bit or 32-bit filterbank files.

All code was developed for e-MERLIN's LOFT-e observing system.

Code was developed by Charles Walker (walker.mbcxqcw2@gmail.com). Median clipping algorithm was developed by Rene Breton.

Please cite accordingly. **(Note to Charlie: Add to ASCL!)**

The RFI clipping algorithm is described in detail in: Charles Walker's thesis: https://www.research.manchester.ac.uk/portal/en/theses/localising-fast-transients(d2847c6c-1b0f-43d9-9155-f8dccb432cb9).html

**(Note to charlie: add this to arXiv!)**

The data file "Cm1670.File1.8bit.fil" referred to in this README is an observation of PSRB0329+54 taken using the LOFT-e backend. It can be found at: https://drive.google.com/open?id=195SnkU3mZbj2cZrPNAwkP5AfGToLVzxT

The ephemeris used for analysing this pulsar is taken from the ATNF pulsar catalog
Manchester, R. N., Hobbs, G.B., Teoh, A. & Hobbs, M., AJ, 129, 1993-2006 (2005)
http://www.atnf.csiro.au/research/pulsar/psrcat 


# PYTHON REQUIREMENTS

os
numpy
matplotlib
astropy
sigpyproc

# EXTRA REQUIREMENTS

dspsr (for folding pulsar data)
pav   (for plotting archive files)

# FILE REQUIREMENTS

"Cm1670.File1.8bit.fil" (a filterbank file to be used to test the code)
"Mary_test.fil"         (the result of running code on "Cm1670.File1.8bit.fil")

# INSTRUCTIONS:

1) Download "Cm1670.File1.8bit.fil" and place it in this directory
2) Download "Mary_test.fil" and place it somewhere else for reference

3) Run "python Mary_Clipping_Script.py" from the command line
4) The file "Cm1670.File1.8bit.fil" will be RFI clipped and output as "Mary_test.fil"
5) Compare the output file to the reference file downloaded. These should be identical.

6) Run dspsr commands on the command line to generate plots to compare PSRB0329+54 pulse profile pre- and post- RFI clipping. This will allow you to compare the two filterbank files and see if the RFI mitigation worked:

6a)   dspsr Cm1670.File1.8bit.fil -F 256 -L 10 -t 8 -k 'Jodrell' -b 1024 -E 0329+54.eph -A -O Mary_RFI

6b)   dspsr Mary_test.fil -F 256 -L 10 -t 8 -k 'Jodrell' -b 1024 -E 0329+54.eph -A -O Mary_clean

6c)   pav -N 1,2 -DFTp Mary*.ar

