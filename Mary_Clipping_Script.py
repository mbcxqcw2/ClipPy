#imports
from clip_filterbanks_2 import ClipFil
from clip_filterbanks_2 import ClipFilFast

#data location and file name
filfile_toclip = 'Cm1670.File1.8bit.fil'

#output clipped file details
out_path = '.'
outfile_name = 'Mary_test.fil'
outfile_fast_name = 'Mary_fast.fil'


#clipping options
bitswap      = False #we wish to keep the output data bit rate the same as the input data bit rate
sigclip      = 3.    #the clipping threshold to use when clipping RFI
toload_samps = 40000 # the number of samples to load at once while clipping

#Perform RFI clipping
ClipFil(filfile_toclip,outfile_name,out_path+'/',bitswap,True,sigclip,toload_samps)
ClipFilFast(filfile_toclip,outfile_fast_name,out_path+'/',bitswap,True,sigclip,toload_samps)
