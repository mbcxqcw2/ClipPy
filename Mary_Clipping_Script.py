#imports
from clip_filterbanks_2 import ClipFil

#data location and file name
filfile_toclip = 'Cm1670.File1.8bit.fil'

#output clipped file details
out_path = '.'
outfile_name = 'Mary_test.fil'


#clipping options
bitswap = False #we wish to keep the output data bit rate the same as the input data bit rate
sigclip = 3.    #the clipping threshold to use when clipping RFI

#Perform RFI clipping
ClipFil(filfile_toclip,outfile_name,out_path+'/',bitswap,True,sigclip)

#running the following dspsr commands will allow you to compare the two filterbank files to see if the RFI mitigation worked:

#   dspsr Cm1670.File1.8bit.fil -F 256 -L 10 -t 8 -k 'Jodrell' -b 1024 =E 0329+54.eph -O Mary_RFI
#   dspsr Mary_test.fil -F 256 -L 10 -t 8 -k 'Jodrell' -b 1024 -E 0329+54.eph -O Mary_clean
#   pav -N 1,2 -DFTp Mary*.ar
