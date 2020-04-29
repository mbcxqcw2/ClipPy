# PREAMBLE

This script performs time-domain RFI clipping on 8-bit or 32-bit filterbank files.

RFI clipping is run twice: once with the single-CPU ClipFil() function, and then again with its multiprocessing ClipFilFast() version.

All code was developed for e-MERLIN's LOFT-e observing system.

Code was developed by Charles Walker (walker.mbcxqcw2@gmail.com). Median clipping algorithm was developed by Rene Breton. Code to initialise a new filterbank file uses the read_header() function developed by Scott Ransom for PRESTO: https://github.com/scottransom/presto/blob/4959ffeac8d369d8883987287179297f93aea773/python/presto/filterbank.py

The RFI clipping algorithm is described in detail in: Charles Walker's thesis: https://www.research.manchester.ac.uk/portal/en/theses/localising-fast-transients(d2847c6c-1b0f-43d9-9155-f8dccb432cb9).html

**(Note to charlie: add this to arXiv!)**

The data file "Cm1670.File1.8bit.fil" referred to in this README is an observation of PSRB0329+54 taken using the LOFT-e backend. It can be found at: https://drive.google.com/open?id=195SnkU3mZbj2cZrPNAwkP5AfGToLVzxT

The ephemeris used for analysing this pulsar is taken from the ATNF pulsar catalog
Manchester, R. N., Hobbs, G.B., Teoh, A. & Hobbs, M., AJ, 129, 1993-2006 (2005)
http://www.atnf.csiro.au/research/pulsar/psrcat 

Please cite accordingly. **(Note to Charlie: Add to ASCL?)**


# PYTHON REQUIREMENTS

os

numpy

matplotlib

multiprocessing

contextlib

astropy

sigpyproc

presto (version 3.0+)

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
5) The file will also be RFI clipped and output by the multiprocessor version of the code and output as "Mary_fast.fil".
6) Compare the output files to the reference file downloaded. These should be identical.

7) Run dspsr commands on the command line to generate plots to compare PSRB0329+54 pulse profile pre- and post- RFI clipping. This will allow you to compare the three filterbank files and see if the RFI mitigation worked:

7a)   dspsr Cm1670.File1.8bit.fil -F 256 -L 10 -t 8 -k 'Jodrell' -b 1024 -E 0329+54.eph -A -O Mary_RFI

7b)   dspsr Mary_test.fil -F 256 -L 10 -t 8 -k 'Jodrell' -b 1024 -E 0329+54.eph -A -O Mary_clean

7c)   dspsr Mary_fast.fil -F 256 -L 10 -t 8 -k 'Jodrell' -b 1024 -E 0329+54.eph -A -O Mary_fastclean

7d)   pav -N 1,3 -DFTp Mary*.ar --ch 2

7e)   pav -N 1,3 -YFd Mary*.ar --ch 2

7f)   pav -N 1,3 -GTd Mary*.ar --ch 2

