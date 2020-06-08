#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
clip_filterbanks_2.py

code to clip filterbank files based on data location directory, mask location, and file extension.
Based on incoherent beam functions in Incoherent_7.py (/share/nas1/LOFT-e/software/incoherent_beam_pipeline).

@author: cwalker (walker.mbcxqcw2@gmail.com)

V1: 20180802

V2: 20180803  - Amended following experiments in:
               /share/nas1/LOFT-e/experiments/02-08-18_rfi_mitigation_test_BWS.
              - Updated Median_clip() to use np.ma.median, allowing the median of
                masked arrays to be calculated correctly
              - Updated RFIclip() to not divide channels by zero if the standard
                deviation of the channel is zero (which follows through with nan errors)
              - Commented out warning when replacing a saturated channel (10/09/2018).
                this is because when running on filmerged filterbanks you get loads of
                printouts. (lines 467 and 506)

V3: 20191219  - Fixed ClipFil() description.
              - Made toload_samps an input variable in ClipFil().
              - Fixed hardcoded RFIClip()instances of drawing numbers
                from Gaussians with means of of X/256 to X/np.float(nchans).
              - Functionised rescaling portion of RFIClip() as RescaleChunk().
              - Functionised cleaning portion of RFIClip() as CleanChunk().
              - Removed superflous loop over # telescopes in ClipFil().
              - Removed dependance of ClipFil() on superflous Beam() function.
              - Removed superflous Beam() function.
              - Removed dependance on superflous RFIClip() function.
              - Removed superflous RFIClip() function.
              - Tidied up structure and comments.
              - Added Multiprocessing ClipFilFast() function.

V4: 20200106  - Tidied up ClipFilFast() function.
              - Made ClipFilFast() dependent on rficlip option again.
              - Cleaned up ReadChunk() helper function.
              - Cleaned up RecastChunk() helper function.
              - Cleaned up RescaleChunk_unwrap() helper function.
              - Cleaned up CleanChunk_unwrap() helper function.
              - Cleaned up RecastChunk_unwrap() helper function.
              - Removed uneccesary import of itertools.product
              - Amended print statements.

V5: 20200107  - Amended ClipFil() to process remainder timesamples rather
                than skipping them. Note: this has since been commented out, 
                because while it works correctly, it rather drastically lowers
                the signal to noise of the pulsar in the test observation.
                More tests on more RFI are neeeded to determine whether it is
                worth keeping the remaining portion. It may be that proper
                statistics cannot be calculated on it as it is so small.
    20200108  - Put ClipFil() remainder timesample code in optional loop
                triggered by "proc_remainder" variable to allow testing.
              - Added proc_remainder option to ClipFilFast() for testing.
    20200211  - Stopped ClipFilFast() from defaulting to using all cores on
                machine

V6: 20200429  - Amended print statements to follow Python 3 conventions.
              - Amended presto sigproc import method to match most recent 
                (Python-3 compatable) version.
              - Modified  ClipFil() to output a standard deviation log file.
              - Modified RescaleChunk() to output a list of standard deviations
                calculated for the channels in its input data prior to rescaling.
              - Amended instances calling RescaleChunk() to handle extra output.
              - Modified ClipFil() to write standard deviations to log file.
    20200806  - Modified Median_clip() for Python 3 compatability using 
                xrange->range conversion where required
              - Amended ClipFilFast() n_complete_chunk_sets calculation to explicitly
                declare integer division (//) to ensure Python 3 compatability

              
"""

from sigpyproc.Readers import FilReader as fr
from astropy.time import Time
from astropy import units as u
import numpy as np
from presto import sigproc
import sigpyproc.Header as spph
import sigpyproc.Utils as sppu
from matplotlib import pyplot as plt
from numpy.random import normal as norm
import os
import multiprocessing as m
from multiprocessing import Pool
from contextlib import closing

def CombineFilUtils_FBchunking(outsamps,blocksize=1000):
    """
    Calculates the number of chunks the filterbank will be read out in
    based on the total number of output samples and the desired
    block reading size.
    
    INPUTS:
    
    outsamps  : (int) the number of samples in the filterbank being written
    blocksize : (int) the number of samples to read at a time
    
    RETURNS:
    
    nchunks   : (int) the number of chunks the filterbank will be written in
    remainder : (int) the remaining number of timesamples after writing final
                      whole block
    
    """
    
    nchunks=int(outsamps/blocksize) #whole number of blocks read in total
    remainder=int(outsamps)%blocksize #remainder of timesamples after final whole block
    
    return nchunks,remainder

def CombineFilUtils_InitialiseOut(fil_names,outloc,outname,startTime,bitswap=False):
    """
    Initialises the output combined filterbank file.
    
    INPUTS:
    
    fil_names : (list) names of filterbanks to combine
    outloc    : (str) desired output filterbank location
    outname   : (str) desired output filterbank file name
    startTime : (mjd) the start time of the output filterbank
    
    RETURNS:
    
    fh_out  : output filterbank file handle
    bitrate : output filterbank bit bits per sample
    
    """

    #get template header information from first filterbank
    head_dict,headsize = read_header(fil_names[0])

    #set output path
    fln_out = outloc + outname
    
    #initialise new header
    header={}
    
    #create and update header
    if bitswap==False:
        fb_hdr = Make_header(sname='junk:incoherent',
                             nbits=head_dict['nbits'],
                             tsamp=head_dict['tsamp'],
                             nchans=head_dict['nchans'],
                             tstart=startTime.mjd,
                             nifs=head_dict['nifs'],
                             fch1=head_dict['fch1'],
                             foff=head_dict['foff'],  **header)
        
        #grab filterbank bit size
        bitrate = head_dict["nbits"]
        
    elif bitswap==True:
        
        #swap the filterbank header bit size from 8 to 32 or vice versa
        filbit=head_dict['nbits']
        if filbit==32:
            filbit=8
        elif filbit==8:
            filbit=32
            
        fb_hdr = Make_header(sname='junk:incoherent',
                             nbits=filbit,
                             tsamp=head_dict['tsamp'],
                             nchans=head_dict['nchans'],
                             tstart=startTime.mjd,
                             nifs=head_dict['nifs'],
                             fch1=head_dict['fch1'],
                             foff=head_dict['foff'],  **header)
        
        #grab filterbank bit size
        bitrate = filbit
    
    #initialise output file handle
    fh_out = []
    
    #initialise filterbank with appropriate header info
    fh_out.append( Init_filterbank(fln_out, fb_hdr) )
    
    return fh_out,bitrate


def Make_header(nbits=32,
                tsamp=1,
                nchans=256,
                tstart=0,
                tel_id=82,
                mac_id=82,
                d_type=1,
                raw="foo.txt",
                sname="junkpsr",
                bcent=0,
                pcent=0,
                az_s=0,
                za_s=0,
                raj=0,
                dej=0,
                ra='00:00:00.000',
                dec='00:00:00.000',
                nsamp=1,
                fch1=1000,
                foff=1,
                nifs=1):
    """
    Make a sigpyproc header based on input header elements.
    Python floats are the equivalent of c doubles, so floats **should** work.
    """
    ## make header data dictionary
    tmp_hdr = {
        "telescope_id":int(tel_id),
        "machine_id":int(mac_id),
        "data_type":int(d_type),
        "rawdatafile":str(raw),
        "source_name":str(sname),
        "barycentric":int(bcent),
        "pulsarcentric":int(pcent),
        "az_start":float(az_s),
        "za_start":float(za_s),
        "src_raj":float(raj),
        "src_dej":float(dej),
        "ra":ra,
        "dec":dec,
        "tstart":float(tstart),
        "tsamp":float(tsamp),
        "nbits":int(nbits),
        "nsamples":int(nsamp),
        "fch1":float(fch1),
        "foff":float(foff),
        "nchans":int(nchans),
        "nifs":int(nifs)
    }

    ## make the sigpyproc header
    hdr = spph.Header(tmp_hdr)

    return hdr

def Init_filterbank(fln, hdr):
    """
    Make a filterbank file and initialise it with the provided header keywords.

    Parameters
    ----------
    fln : str
        Filename of the filterbank.
    **kwargs
        Keywords inputs from the Make_header function.

    Return
    ------
    sigpyproc filterbank file handler.
    """
    fh = spph.Header.prepOutfile(hdr, fln)

    return fh    

def read_header(filename, verbose=False):
    """Read the header of a filterbank file, and return
        a dictionary of header paramters and the header's
        size in bytes.
        Inputs:
            filename: Name of the filterbank file.
            verbose: If True, be verbose. (Default: be quiet)
        Outputs:
            header: A dictionary of header paramters.
            header_size: The size of the header in bytes.

    Note: this was borrowed from Scott Ransom's presto github:

    https://github.com/scottransom/presto/blob/master/lib/python/filterbank.py

    Note Note: this is a direct copy of the read_header in downsamp_utils.
    I should import it really...

    """
    header = {}
    print('to read: '+filename)
    filfile = open(filename, 'rb')
    filfile.seek(0)
    paramname = ""
    while (paramname != 'HEADER_END'):
        if verbose:
            print("File location: %d" % filfile.tell())
        paramname, val = sigproc.read_hdr_val(filfile, stdout=verbose)
        if verbose:
            print("Read param %s (value: %s)" % (paramname, val))
        if paramname not in ["HEADER_START", "HEADER_END"]:
            header[paramname] = val
    header_size = filfile.tell()
    filfile.close()
    return header, header_size



def DownSampleBits(data,clip=4):
    """
    
    Downsamples 2d array (i.e. a filterbank file) of 32-bit floats.
    Returns an 8-bit array.
    
    INPUT:
    
    data : (array-like) The array to downsample
    clip : (integer)    Clipping value. Anything this number
                        of standard deviations away from
                        the mean will be clipped.
                        
    RETURNS:
    
    data : (array-like) 8-bit downsampled data.
    
    """
    
    data-=np.mean(data)
    #print 'ERROR CHECK: ', data,np.mean(data),np.std(data),' END ERROR CHECK'
    data/=np.std(data)
    data*=128./clip
    data+=128
    data=data.astype(np.uint8)
    
    return data


def CombineFilUtils_FBoverlap(fils):
    """
    Calculates how much overlap in time exists between a set of
    filterbank files to combine, and what time samples
    to skip when combining them.
    
    INPUTS:
    
    fils : (list) list of pointers to the filterbank files
                  to combine, as read in by sigpyproc.FilReader
                  
    RETURNS:
    
    outsamps : (int) the number of samples that will be read
                     for each filterbank
                     
    nskips   : (list) list containing the number of samples to
               skip from the beginning of each input filterbank
               when combining them.
               
    t_i      : (mjd) largest start time from the list of filterbanks,
                     which will be the start time of the combined
                     filterbank.
                     
    nchans   : (int) the number of channels in the filterbanks to combine
    
    """
    
    #number of filterbanks
    nFils = len(fils)
    
    #initialise filterbank info arrays
    
    nchans = []  #number of channels
    tstarts = [] #start times of filterbanks
    nsamps = []  #number of samples in filterbanks
    tsamps = []  #sampling times of filterbanks
    dtMax = []   #maximum time offset from tstart for each filterbank
    
    #fill arrays with appropriate filterbank information
    
    for fil in fils:                                           #loop over filterbanks
        
        nchans.append(fil.header.nchans)                       #number of channels
        tstarts.append(Time(fil.header.tstart,format='mjd'))   #start mjd
        nsamps.append(fil.header.nsamples)                     #number of time samps
        tsamps.append(fil.header.tsamp)                        #sampling time [seconds]
        dtMax.append(u.s*fil.header.nsamples*fil.header.tsamp) #max offset from tstart [s]
        
    #get all filterbank end times
    
    tends = [tstarts[i]+dtMax[i] for i in range(nFils)]
    
    #find extremes of filterbanks for cropping purposes
    
    i = np.argmax(tstarts) #index of filterbank with largest start time
    t_i = tstarts[i]       #largest start time
    e = np.argmin(tends)   #index of filterbank with smallest end time
    t_e = tends[e]         #smallest end time
    
    #find total number of samples which will be read based on extremes
    
    outsamps = (t_e.tt - t_i.tt).to('s')/(tsamps[0]*u.s) #output fb nsamps = time difference / sampling time

    #calculate the number of samples to skip at beginning of each filterbank
    
    nskips = [(t_i.tt-tstarts[i].tt).to('s')/(tsamps[0]*u.s) for i in range(len(tstarts))] #number of time samples to skip from each file

    return outsamps,nskips,t_i,nchans[0]

###########################################################
##### MEDIAN CLIPPING FUNCTIONS FOR RFI MITIGATION ########
###########################################################

def Median_clip(arr, sigma=3, max_iter=3, ftol=0.01, xtol=0.05, full_output=False, axis=None):
    """

    #from /home/bretonr/lib/python/pyastro/misc.py

    #modified by cwalker on 08/06/2020 to support Python 3 compatability

    Median_clip(arr, sigma, max_iter=3, ftol=0.01, xtol=0.05, full_output=False, axis=None)
    Return the median of an array after iteratively clipping the outliers.
    The median is calculated upon discarding elements that deviate more than
    sigma * standard deviation the median.
    
    arr (numpy.array): array to calculate the median from.
    sigma (float): the clipping threshold, in units of standard deviation.
    max_iter (int): the maximum number of iterations. A value of 0 will
        return the usual median.
    ftol (float): fraction tolerance limit for convergence. If the number
        of discarded elements changes by less than ftol, the iteration is
        stopped.
    xtol (float): absolute tolerance limit for convergence. If the number
        of discarded elements increases above xtol with respect to the
        initial number of elements, the iteration is stopped.
    full_output (bool): If True, will also return the indices that were good.
    axis (None/int): Axis along which the calculation is to be done. NOT WORKING!!!
    
    >>> med = Median_clip(arr, sigma=3, max_iter=3)
    >>> med, std, inds_good = Median_clip(arr, sigma=3, max_iter=3, full_output=True)
    
    med (float): clipped median.
    std (float): clipped standard deviation.
    inds_good (numpy.array): array of bool where True values are those used to
        compute the clipped median.
    """
    arr = np.ma.masked_invalid(arr)
    med = np.ma.median(arr, axis=axis)
    std = np.ma.std(arr, axis=axis)
    ncount = arr.count(axis=axis)

    #perform a check to see if python 2.7 or free is running, and use the appropriate (x)range function
    try:
        #Python 2
        xrange
    except NameError:
        #Python 3, xrange is now named range
        xrange=range


    for niter in xrange(max_iter):
        ncount_old = arr.count(axis=axis)
        if axis is not None:
            condition = (arr < np.expand_dims(med-std*sigma, axis)) + (arr > np.expand_dims(med+std*sigma, axis))
        else:
            condition = (arr < med-std*sigma) + (arr > med+std*sigma)
        arr = np.ma.masked_where(condition, arr)
        ncount_new = arr.count(axis)
        med = np.ma.median(arr, axis=axis)
        std = np.ma.std(arr, axis=axis)
        #print med,std
        if np.any(ncount-ncount_new > xtol*ncount):
            #print( "xtol reached {}; breaking at iteration {}".format(1-1.*ncount_new/ncount, niter+1) )
            break
        if np.any(ncount_old-ncount_new < ftol*ncount_old):
            #print( "ftol reached {}; breaking at iteration {}".format(1-1.*ncount_new/ncount_old, niter+1) )
            break
        if niter == max_iter-1:
            print( "maximum number of iterations reached" )
    if full_output:
        if isinstance(arr.mask, np.bool_):
            mask = np.ones(arr.shape, dtype=bool)
        else:
            mask = ~arr.mask
        if axis is not None:
            med = med.data
            std = std.data
        return med, std, mask
    if axis is not None:
        med = med.data
    return med

def RescaleChunk(data,nchans,sig=3.):
    """
    Individually rescales each channel in a chunk of data to have mean 0 and std 1.
    Calls iterative median clipping algorithm Median_clip().

 
    INPUTS:

    data : (array-like) filterbank data chunk
    nchans : (int) number of filterbank channels in original file
    sig: (float) standard deviations away from mean to clip after


    RETURNS:

    rescaled_data : (array-like) rescaled data chunk
    std_list : (list) list of the standard deviations calculated for and used to
    rescale each channel in the data chunk by the median clipping algorithm. If
    data was readin using sigpyproc.readBlock() std list will follow the same convention
    lowest frequency channel -> highest frequency channel.

    """

    std_list = [] #initialise list to hold channel standard deviations

    for j in range(nchans):
        #data
        channel = data[j]
        #get mean, standard deviation
        mean,std,mask = Median_clip(channel,sig,max_iter=5,full_output=True,xtol=0.5)#edit: max_iter 10>5
        #record channel standard deviation
        std_list.append(std)
        if std==0.0:
            #print 'WARNING: Channel', j, 'mean, std ',mean,std
            channel-=mean
        else:
            channel-=mean
        #divide by std
            channel/=std

    rescaled_data = data

    return rescaled_data,std_list


def CleanChunk(data,nchans,sig=3.):
    """
    Clips time-domain RFI from a chunk of filterbank data which has been cleaned
    using RescaleChunk().
    Calls iterative median clipping algorithm Median_clip().


    ALGORITHM:

    1) Crunch data to get timeseries

    2) Get median and stdv of timeseries

    3) Find where timeseries lies outside of predefined sigma level

    4) On channel-by-channel basis replace bad timesamples with random numbers drawn from gaussian


    INPUTS:

    data   : (array-like) filterbank data chunk which has been rescaled using RescaleChunk()
    nchans : (int) number of filterbank channels in original file
    sig    : (float) standard deviations away from mean to clip after


    RETURNS:

    cleaned_data : (array-like) RFI-clipped data chunk
    """

    #make timeseries from data chunk
    timeseries=data.sum(axis=0)

    #get median and std of the timeseries
    # the sigma for median clipping is currently the same as the rescaling
    
    med,std,mask = Median_clip(timeseries,sig,max_iter=5,full_output=True,xtol=0.5) #edit: max_iter 10>5
    #print 'timeseries mean, std ', med,std

    #find where timeseries data lies outside of boundaries
    #boundaries are std standard deviaions away from
    #the timeseries median
    
    minval = med-(sig*std)
    maxval = med+(sig*std)
    #print 'minval, maxval ',minval,maxval
    toclip = ((timeseries<minval)|(timeseries>maxval))

    #loop over data on channel by channel basis
    #for each channel, replace timesamples which, in the timeseries, lay outside
    #of the clipping range, with new numbers
    #the numbers are randomly drawn from a gaussian
    #the gaussian has a mean of the timeseries median/nchans (so when summed, they will
    #lie around the correct mean) and a standard deviation of the channel (which ideally
    #should be 1)
    
    for i in range(data.shape[0]):
        #select channel
        channel=data[i,:]
        #get good median and std
        chan_med,chan_std,mask = Median_clip(channel,sig,max_iter=5,full_output=True,xtol=0.5)#edit: max_iter 10>5
        if chan_std<=0:
            #print 'WARNING: Channel', i, 'mean, chan_std ',chan_med,chan_std

            # change: 19/07/2018
            # if this warning occurs, the entire channel is
            # saturated for this timeseries and should be replaced with randoms
            # drawn from gaussian with mean: median/nchans and std: 1 (I think) so insert below line:
            chan_std = 1.
            channel=np.random.normal(loc=med/np.float(nchans),scale=1,size=channel.shape)

        #replace bad data with median of timeseries/nchans and std of channel std
        channel = (toclip*np.random.normal(loc=med/np.float(nchans),scale=chan_std,size=channel.shape))+(~toclip*channel)
        #overwrite
        data[i,:]=channel

    #rename cleaned array
    cleaned_data = data

    return cleaned_data


#################################################################################

def ClipFil(in_fil,outname,outloc,bitswap,rficlip=True,clipsig=3.,toload_samps=40000,proc_remainder=False):
    """
    Mitigates timeseries RFI in an input filterbank file.
    Outputs second, clipped file with chosen name.
    Mitigates RFI on a chunk-by-chunk basis using a chosen chunksize.

    Also outputs the log file "stdlog.txt" which returns the standard
    deviations calculated for each frequency channel and data chunk
    when rescaling original data. Structure of stdlog.txt is:

    Line 1: header
    Line 2: chunk 1: standard deviation (lowest frequency) -> standard deviation (highest frequency)
    Line 3: chunk 2: standard deviation (lowest frequency) -> standard deviation (highest frequency)
    ...
    etc to final chunk.


    CHUNKWISE RFI MITIGATION ALGORITHM:

    1) Reads chunk

    2) Calls RescaleChunk()

        a) Individually rescales channels to have mean 0 and stdv 1.

    3) Calls CleanChunk()

        b) Crunch rescaled data to get timeseries

        c) Get median and stdv of timeseries

        d) Find where timeseries lies outside of predefined sigma level

        e) On channel-by-channel basis replace bad timesamples with random numbers drawn from gaussian

    4) Returns rescaled, cleaned data chunk


    FUNCTION INPUTS:

    in_fil         : (str) input filterbank to clip (must be directory location and
                     file name)
    outname        : (str) output name for clipped filterbank
    outloc         : (str) output folder for clipped filterbank
    bitswap        : (True/False boolean) if True, 8-bit input will be written out
                     as 32-bit
                     if False, 8-bit input will be written out as 8-bit
                     and vice-versa
    rficlip        : (True/False boolean) if True, rfi sigma clipping will be applied via
                     timeseries data
    clipsig        : (float) the clipping sigma, if rficlip==True
    toload_samps   : (int) the number of timesamples (spectra) to load at once while
                     clipping.
    proc_remainder : (boolean) option to process any remaining timesamples which do not
                     fit into an integer number of toload_samps. If "False", these
                     timesamples are discarded. Default is "False" until
                     testing of its effect on S/N is complete.



    """
    

    ##INITIALISE INPUT FILTERBANK LOADING##
    print('Input file is: ',in_fil)
    fil_names=np.array([in_fil]).flatten()
    print(fil_names.shape)
    print('loading filterbank\n')
    fils=[]
    for fil in fil_names:
        fils.append(fr(fil)) #store pointers to filterbank file
    print('Calculating start mjd and samples to read and skip\n')
    #note, nskips should be [0] and outsamps should be the number of
    # timesamples in the filterbank file as there is only one input.

    outsamps,nskips,startTime,nchans = CombineFilUtils_FBoverlap(fils)
    nskips=np.zeros_like(nskips)
    print('    ...samples to skip (should be [0]): {0}'.format(nskips))
    print('    ...samples to read (should be fil length): {0}'.format(outsamps))



    ##INITIALISE CHUNKING INFORMATION##
    print('calculating data chunking information\n')
    blocksize=toload_samps
    nchunks,remainder = CombineFilUtils_FBchunking(outsamps,blocksize)
    


    ##INITIALISE THE OUTPUT FILTERBANK FILE
    print('initialising output filterbank\n')
    fh_out,bitrate = CombineFilUtils_InitialiseOut(fil_names,
                                                   outloc,
                                                   outname,
                                                   startTime,
                                                   bitswap=bitswap)
    


    ##SET DATA TYPE TO WRITE OUTPUT AS##
    print('output files will be {0}-bit\n'.format(bitrate))
    if bitrate==8:
        outdtype = np.uint8
    elif bitrate==32:
        outdtype = np.float32



    ##PERFORM RFI MITIGATION AND WRITE TO FILE##

    #open standard deviation log file and write header
    stdlog = open("{0}/stdlog.txt".format(outloc),"w")
    stdlog.write("Original number of time samples: {0} Number of frequency channels: {1} Processing chunk size (timesamps): {2} Number of chunks: {3} Process remainder: {4} Clipping sigma: {5} File structure: rows=chunks, columns=frequency channel standard deviations (low freq->high freq)\n".format(outsamps,nchans,toload_samps,nchunks,proc_remainder,clipsig))

    #begin clipping
    print('beginning clipping\n')
   
    nfils = len(fils) #number of filterbank files to clip. Should always be 1

    for c in range(nchunks): #loop over chunks to write out

        ###INITIALISE CHUNK LOADING###
        data=np.zeros((nchans,blocksize,nfils)) #declare 3D array to hold data
        chunk = 0 #initialise filterbank chunk
        skip=int(round(nskips[0]))
        #Note: number of blocks to skip reading at beginning of file (nskips[0] should
        # be zero, this is vestigial from incoherent beam code

        blockstart=int(skip+(c*blocksize)) #start sample of chunk to read

        print('Reading/Writing chunk {0}/{1}'.format(c,nchunks))

        ###READ CHUNK###
        chunk=fils[0].readBlock(blockstart,blocksize) #read chunk

        ###OPTIONAL: RESCALING AND CLIPPING###
        if rficlip==True: #if rfi clipping mode is on:
            print('RFI clipping...')

            ###RESCALE CHUNK###
            chunk,stdlist=RescaleChunk(chunk,nchans,clipsig) #rescale the chunk
            stdlog.write(" ".join(np.array(stdlist,dtype=str))+"\n") #write standard deviations to log file

            ###CLIP CHUNK###
            chunk=CleanChunk(chunk,nchans,clipsig)

        ###STORE CLEANED, RESCALED CHUNK IN NEW ARRAY###
        data[:,:,0]=chunk #append telescope to data
        data=data.sum(axis=2)
        #flatten data array over third axis. Transforms data from 3D array
        # (of shape: [channels,times,1]) to 2D array (of shape: [channels,times])

        ###OPTIONAL: RESCALING OF DATA PRODUCT FOR STORAGE###
        if bitrate==8: #if necessary...
            data=DownSampleBits(data) #...downsample to 8-bit

        ###WRITE OUT CLEANED DATA TO FILE###
        data=data.T.flatten().astype(dtype=outdtype)
        #reshape the data to filterbank output (low freq to high freq t1, low
        # freq to high freq t2, ....) and recast to desired type

        sppu.File.cwrite(fh_out[0], data) #write block to filterbank file


    ##PROCESS REMAINING SAMPLES##
    ### NOTE: this defaults to false pending further investigation (see version notes V5) ###
    if proc_remainder==True:
        print('Remaining samples will be clipped. Clipping remaining {0} samples...'.format(remainder))
    
        #initialise chunk load
        data = np.zeros((nchans,remainder,nfils))
        skip = int(round(nskips[0])) #should always be zero
        blockstart = int(skip+outsamps-remainder+1)
        #NOTE: skip should always be zero. This is a vestigial variable from
        # incoherent beamforming code. Outsamps is the number of time samples
        # in the filterbank file. Python is zero-indexed, hence you must
        # start loading from the sample (outsamps - remainder + 1).
    
        #read remainder
        print('    Reading remainder...')
        chunk = fils[0].readBlock(blockstart,remainder)
    
        #optional: rescaling and clipping
        if rficlip==True: #if rfi clipping mode is on:
            print('    Rescaling remainder...')
            chunk,stdlist = RescaleChunk(chunk,nchans,clipsig) #rescale the chunk
            stdlog.write(" ".join(np.array(stdlist,dtype=str))+"\n") #write standard deviations to log file
            print('    Cleaning remainder...')
            chunk = CleanChunk(chunk,nchans,clipsig)
    
        #store clean, rescaled chunk in new array
        data[:,:,0]=chunk
        data=data.sum(axis=2)
    
        #optional: rescale data product for storage
        if bitrate==8:
            print('    Downsampling remainder...')
            data=DownSampleBits(data)
    
        #write clean data to file
        data = data.T.flatten().astype(dtype=outdtype)
        #reshape the data to filterbank output (low freq to high freq t1,
        # low freq to high freq t2, ....) and recast to desired type
    
        sppu.File.cwrite(fh_out[0],data) #write remainder to filterbank file


    stdlog.close()
    print('ClipFil() process complete.')
    ##END FUNCTION##



    return



def ReadChunk(datachunk):
    """
    A helper function for reading data chunk via multiprocessing. Exists because you
    can't multiprocess using a function with multiple inputs in Python 2.7 very easily.
    See stackoverflow thread: "Python multiprocessing pool.map for multiple arguments".


    INPUTS:

    datachunk : (array-like) a chunk of filterbank data read by sigpyproc's readBlock()

    RETURNS:

    data      : (array-like) a copy of datachunk

    """

    data=datachunk#fils[0].readBlock(startsample,blocksize)

    return data

def RecastChunk(datachunk,outdtype):
    """
    A helper function for recasting data chunk via multiprocessing. Exists because you
    can't multiprocess using a function with multiple inputs in Python 2.7 very easily.
    See stackoverflow thread: "Python multiprocessing pool.map for multiple arguments".

    This function reshapes input data to sigpyproc filterbank output (low freq to
    high freq t1, low freq to high freq t2, ....) and recasts it to a desired
    output bit float type.

    INPUTS:

    datachunk : (array-like) a chunk of filterbank data read by sigpyproc's readBlock().
    outdtype  : either np.uint8 or np.float32

    RETURNS:

    datachunk : (array-like) recast/reshaped filterbank data.

    """

    datachunk=datachunk.T.flatten().astype(dtype=outdtype)

    return datachunk

def RescaleChunk_unwrap(args):
    """
    A helper function for rescaling data chunk via multiprocessing. Exists because
    you can't multiprocess using a function with multiple inputs in Python 2.7 very
    easily. See stackoverflow thread: "Python multiprocessing pool.map for multiple
    arguments".

    This function unpacks a chunk of data and information necessary for rescaling, and
    passes them into RescaleChunk().

    INPUTS:

    args : a list containing 3 arguments in the following order:

           datachunk : (array-like) chunk of filterbank data read by sigpyproc's readBlock()
           nchans    : (int) number of filterbank channels in original file
           sig       : (float) standard deviations away from mean to clip after

    RETURNS:

    output of RescaleChunk() (Currently only the data and not the list of standard deviations)

    """

    datachunk=args[0]
    nchans=args[1]
    sig=args[2]

    return RescaleChunk(datachunk,nchans,sig)[0]

def CleanChunk_unwrap(args):
    """
    A helper function for cleaning data chunk via multiprocessing. Exists because you
    can't multiprocess using a function with multiple inputs in Python 2.7 very easily.
    See stackoverflow thread: "Python multiprocessing pool.map for multiple arguments".

    This function unpacks a chunk of data and information necessary for cleaning, and
    passes them into CleanChunk().

    INPUTS:

    args : a list containing 3 arguments in the following order:

           rescaledchunk : (array-like) filterbank data chunk read in by sigpyproc's
                           readBlock() and rescaled using RescaleChunk().
           nchans        : (int) number of filterbank channels in original file.
           sig           : (float) standard deviations away from mean to clip after.

    RETURNS:

    output of CleanChunk()


    """

    #print args
    rescaledchunk=args[0]
    nchans=args[1]
    sig=args[2]

    return CleanChunk(rescaledchunk,nchans,sig)

def RecastChunk_unwrap(args):
    """
    A helper function for recasting data chunk via multiprocessing. Exists because you
    can't multiprocess using a function with multiple inputs in Python 2.7 very easily.
    See stackoverflow thread: "Python multiprocessing pool.map for multiple arguments".

    This function unpacks a chunk of data and information necessary for recasting, and
    passes them into RecastChunk().

    INPUTS:

    args : a list containing 2 arguments in the following order:

           datachunk : (array-like) a chunk of filterbank data read by
           sigpyproc's readBlock()

           outdtype  : either np.uint8 or np.float32

    RETURNS:

    output of RecastChunk()

    """

    datachunk=args[0]
    outdtype=args[1]

    return RecastChunk(datachunk,outdtype)



def ClipFilFast(in_fil,outname,outloc,bitswap,rficlip=True,clipsig=3.,toload_samps=40000,n_cores=2,proc_remainder=False):
    """
    Same as ClipFil but parallelised for speed. See ClipFil().

    Extra inputs:

    n_cores : [int] number of cores to use in parallel.
                    default is 2. If changed to a number greater than the total
                    number of available cores on the machine, the number used 
                    will default to the maximum number of cores on the machine.
                    Using all cores may be unadvisable, so decide this value wisely.


    """
    

    ##INITIALISE INPUT FILTERBANK LOADING##
    print('Input file is: ',in_fil)
    fil_names=np.array([in_fil]).flatten()
    print(fil_names.shape)
    print('loading filterbank\n')
    fils=[]
    for fil in fil_names:
        fils.append(fr(fil)) #store pointers to filterbank file
    print('Calculating start mjd and samples to read and skip\n')
    #note, nskips should be [0] and outsamps should be the number of 
    # timesamples in the filterbank file as there is only one input.
    outsamps,nskips,startTime,nchans = CombineFilUtils_FBoverlap(fils)
    nskips=np.zeros_like(nskips)
    print('    ...samples to skip (should be [0]): {0}'.format(nskips))
    print('    ...samples to read (should be fil length): {0}'.format(outsamps))



    ##INITIALISE CHUNKING INFORMATION##
    print('calculating data chunking information\n')
    blocksize=toload_samps
    nchunks,remainder = CombineFilUtils_FBchunking(outsamps,blocksize)
    


    ##INITIALISE THE OUTPUT FILTERBANK FILE
    print('initialising output filterbank\n')
    fh_out,bitrate = CombineFilUtils_InitialiseOut(fil_names,
                                                   outloc,
                                                   outname,
                                                   startTime,
                                                   bitswap=bitswap)
    


    ##SET DATA TYPE TO WRITE OUTPUT AS##
    print('output files will be {0}-bit\n'.format(bitrate))
    if bitrate==8:
        outdtype = np.uint8
    elif bitrate==32:
        outdtype = np.float32



    #INITIALISE PARALELLISATION##
    ncpus = m.cpu_count() #count available cpus
    print('Maximum number of cpus which may be used at once is: {0}'.format(ncpus))

    #decide number of cores to use
    if n_cores > ncpus: #case: the number of requested cores is greater than the numbeer of available cores
        print("WARNING! Number of requested cores ({0}) is too large! Will use maximum available.".format(n_cores))
        ncpus = ncpus
    elif n_cores <= ncpus: # case: the number of requested cores is less than/equal to the number of available cores
        ncpus = n_cores 

    print('The number of cpus which will be used at once is: {0}'.format(ncpus))
    print('The total number of subchunks to process is: {0}'.format(nchunks))
    print('The number of subchunks to be processed in one chunk set is: {0}'.format(ncpus))

    n_complete_chunk_sets = nchunks//ncpus #calculate the number of sets of chunks where all cpus will be used to be processed
    print('The number of complete chunk sets to be processed is: {0}'.format(n_complete_chunk_sets))

    n_partial_chunk_set = nchunks%ncpus #calculate the number of cpus which must be used to process any remaining chunks
    print('The number of remainder subchunks to be processed is: {0}\n'.format(n_partial_chunk_set))

    subchunklist = np.arange(nchunks) #list of subchunks to process
    print('The subchunks to be processed are: {0}'.format(subchunklist))

    skip=int(round(nskips[0])) #number of blocks to skip reading at beginning of file (=0)
    blockstartlist = [int(skip+(c*blocksize)) for c in subchunklist] #list of chunk start samples to read
    print('The start samples to be read for each subchunk are: {0}\n'.format(blockstartlist))

    #print [blockstartlist[i+(0*ncpus)] for i in range(ncpus)]


    ##PROCESS COMPLETE CHUNK SETS##
    for count in range(n_complete_chunk_sets): #loop over rounds of full cpu usage
        print('Processing complete {0}-subchunk chunk set {1}/{2}...'.format(ncpus,count+1,n_complete_chunk_sets))

        with closing(Pool(ncpus)) as p:
        #invoke multiprocessing (see stackoverflow threads: "Python 3: does
        # Pool keep the original order of data passed to map?" and: "Python
        # Multiprocessing Lib Error (AttributeError: __exit__)"


            ###READ SUBCHUNKS INTO DIFFERENT CPUS###
            subchunks=p.map(ReadChunk, [(fils[0].readBlock(blockstartlist[i+(count*ncpus)],blocksize)) for i in range(ncpus)],chunksize=1)
            #Note: from stackoverflow threads: "Python 3: does Pool keep the original
            # order of data passed to map?" and "multiprocessing pool.map not processing
            # list in order", the use of pool with the variable "chunksize=1" should force
            # subchunks to be processed/output to new array in order if Python 2.7 works
            # in the same way as Python 3. Need to check to see if this is true.


            ###OPTIONAL: RESCALING AND CLIPPING###
            if rficlip==True: #if rfi clipping mode is on:
                print('    RFI clipping = True.')
                #rescale all subchunks in cpus
                print('    Rescaling subchunks...')
                rescaled_subchunks=p.map(RescaleChunk_unwrap,([(subchunk,nchans,clipsig) for subchunk in subchunks]),chunksize=1)
                #clean all subchunks in cpus
                print('    Cleaning subchunks...')
                cleaned_subchunks=p.map(CleanChunk_unwrap,([(r_subchunk,nchans,clipsig) for r_subchunk in rescaled_subchunks]),chunksize=1)
            elif rficlip==False: #else:
                print('    RFI clipping = False. Data will not be clipped.')
                #do not modify the subchunks
                cleaned_subchunks = subchunks


            ###OPTIONAL: RESCALING OF DATA PRODUCT FOR STORAGE###
            if bitrate==8: #if necessary...
                print('    Rescaling data to 8 bit...')
                out_subchunks=p.map(DownSampleBits,[c_subchunk for c_subchunk in cleaned_subchunks],chunksize=1) #...downsample to 8-bit
            else:
                out_subchunks=cleaned_subchunks


            ###RECAST DATA FOR OUTPUT###
            print('    Recasting data for output....')
            recast_subchunks=p.map(RecastChunk_unwrap,([(o_subchunk,outdtype) for o_subchunk in out_subchunks]),chunksize=1)
            #reshape the data to filterbank output (low freq to high freq t1, low
            # freq to high freq t2, ....) and recast to desired bit float type


            p.terminate()


        ###WRITE OUT DATA TO FILE###
        print('    Writing subchunks...')
        for subchunk in recast_subchunks:
            sppu.File.cwrite(fh_out[0], subchunk) #write subchunk to filterbank file


    ##PROCESS PARTIAL REMAINING SUBCHUNKS##
    print('Processing remaining {0} subchunks...'.format(n_partial_chunk_set))


    with closing(Pool(n_partial_chunk_set)) as p: #invoke multiprocessing
       

        ###READ INDIVIDUAL SUBCHUNKS INTO CPUS###
        subchunks=p.map(ReadChunk, [fils[0].readBlock(blockstartlist[i + (ncpus*n_complete_chunk_sets)],blocksize) for i in range(n_partial_chunk_set)],chunksize=1)


        ###OPTIONAL: RESCALING AND CLIPPING###
        if rficlip==True: #if rfi clipping mode is on:
            print('    RFI clipping = True')
            #rescale all subchunks in cpus
            print('    Rescaling subchunks...')
            rescaled_subchunks=p.map(RescaleChunk_unwrap,([(subchunk,nchans,clipsig) for subchunk in subchunks]),chunksize=1)
            #clean all subchunks in cpus
            print('    Cleaning subchunks...')
            cleaned_subchunks=p.map(CleanChunk_unwrap,([(r_subchunk,nchans,clipsig) for r_subchunk in rescaled_subchunks]),chunksize=1)
        elif rficlip==False: #else:
            print('    RFI clipping = False. Data will not be clipped.')
            #do not modify the subchunks
            cleaned_subchunks = subchunks


        ###OPTIONAL: RESCALING OF DATA PRODUCT FOR STORAGE###
        if bitrate==8: #if necessary...
            print('    Rescaling data to 8-bit...')
            out_subchunks=p.map(DownSampleBits,[c_subchunk for c_subchunk in cleaned_subchunks],chunksize=1) #...downsample to 8-bit
        else:
            out_subchunks=cleaned_subchunks


        ###RECAST DATA FOR OUTPUT###
        print('    Recasting data for output...')
        recast_subchunks=p.map(RecastChunk_unwrap,([(o_subchunk,outdtype) for o_subchunk in out_subchunks]),chunksize=1)
        #reshape the data to filterbank output (low freq to high freq t1, low freq to
        # high freq t2, ....) and recast to desired bit float type


        p.terminate()


    ###WRITE OUT DATA TO FILE###
    print('    Writing subchunks...')
    for subchunk in recast_subchunks:
        sppu.File.cwrite(fh_out[0], subchunk) #write block to filterbank file



    ##PROCESS REMAINING SAMPLES##
    ### NOTE: this defaults to false pending further investigation (see version notes V5) ###
    if proc_remainder==True:
        print('Remaining samples will be clipped. Clipping remaining {0} samples...'.format(remainder))
    
        #initialise chunk load
        data = np.zeros((nchans,remainder,len(fils)))
        skip = int(round(nskips[0])) #should always be zero
        blockstart = int(skip+outsamps-remainder+1)
        #NOTE: skip should always be zero. This is a vestigial variable from
        # incoherent beamforming code. Outsamps is the number of time samples
        # in the filterbank file. Python is zero-indexed, hence you must
        # start loading from the sample (outsamps - remainder + 1).
    
        #read remainder
        print('    Reading remainder...')
        remainder_subchunk = fils[0].readBlock(blockstart,remainder)
    
        #optional: rescaling and clipping
        if rficlip==True: #if rfi clipping mode is on:
            print('    Rescaling remainder...')
            remainder_subchunk,stdlist = RescaleChunk(remainder_subchunk,nchans,clipsig)
            print('    Cleaning remainder...')
            remainder_subchunk = CleanChunk(remainder_subchunk,nchans,clipsig)
    
        #store clean, rescaled chunk in new array
        data[:,:,0]=remainder_subchunk
        data=data.sum(axis=2)
    
        #optional: rescale data product for storage
        if bitrate==8:
            print('    Rescaling remainder to 8-bit...')
            data=DownSampleBits(data)
    
        #write clean data to file
        data = data.T.flatten().astype(dtype=outdtype)
        #reshape the data to filterbank output (low freq to high freq t1,
        # low freq to high freq t2, ....) and recast to desired type
    
        sppu.File.cwrite(fh_out[0],data) #write remainder to filterbank file



    print('ClipFilFast() process complete.')
    ##END FUNCTION##



    return
