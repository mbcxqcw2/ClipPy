#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
clip_filterbanks_2.py

code to clip filterbank files based on data location directory, mask location, and file extension. Based on incoherent beam functions
in Incoherent_7.py (/share/nas1/LOFT-e/software/incoherent_beam_pipeline)

@author: cwalker (walker.mbcxqcw2@gmail.com)

V1: 20180802
V2: 20180803  - Amended following experiments in /share/nas1/LOFT-e/experiments/02-08-18_rfi_mitigation_test_BWS.
              - Updated Median_clip() to use np.ma.median, allowing the median of masked arrays to be calculated correctly
              - Updated RFIclip() to not divide channels by zero if the standard deviation of the channel is zero (which follows through with nan errors)
              -commented out warning when replacing a saturated channel (10/09/2018)
              -this is because when running on filmerged filterbanks you get loads of printouts. (lines 467 and 506)
V3: 20191219  - Fixed ClipFil() description.
              - Made toload_samps an input variable in ClipFil().
              - Fixed hardcoded RFIClip()instances of drawing numbers from Gaussians with means of of X/256 to X/np.float(nchans).
              
"""

from sigpyproc.Readers import FilReader as fr
from astropy.time import Time
from astropy import units as u
import numpy as np
import sigproc
import sigpyproc.Header as spph
import sigpyproc.Utils as sppu
from matplotlib import pyplot as plt
from numpy.random import normal as norm
import os

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

    Note: this was borrowed from Scott Ransom's presto github https://github.com/scottransom/presto/blob/master/lib/python/filterbank.py
    Note Note: this is a direct copy of the read_header in downsamp_utils.
    I should import it really...

    """
    header = {}
    print 'to read: '+filename
    filfile = open(filename, 'rb')
    filfile.seek(0)
    paramname = ""
    while (paramname != 'HEADER_END'):
        if verbose:
            print "File location: %d" % filfile.tell()
        paramname, val = sigproc.read_hdr_val(filfile, stdout=verbose)
        if verbose:
            print "Read param %s (value: %s)" % (paramname, val)
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
    print 'ERROR CHECK: ', data,np.mean(data),np.std(data),' END ERROR CHECK'
    data/=np.std(data)
    data*=128./clip
    data+=128
    data=data.astype(np.uint8)
    
    return data


def Beam(data):
    """
    Creates a beam. As the number of telescopes used in ClipFil() is always 1, this is 
    unecessarily complicated and should be cleaned up, as you don't need to crunch over 1 telescope.
    This is basically a vestigal algorithm from Incoherent_7.
    
    INPUT:
    
    data : (array-like) 3D array (of shape: [channels,times,ntelescopes])
                        to beam.
                        
    RETURNS:
    
    ibeam : (array-like) 2D array (of shape: [channels,times])
                         which has been crunched over ntelescopes (which should be 1).
    
    """
    
    ibeam = data.sum(axis=2)
    
    return ibeam

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

def RFIclip(data,nchans,sig=3.):
    """
    Mitigates timeseries RFI in filterbank data.

    ALGORITHM:

    1) Individually rescales channels to have mean 0 and stdv 1

    2) Crunch data to get timeseries

    3) Get median and stdv of timeseries

    4) Find where timeseries lies outside of predefined sigma level

    5) On channel-by-channel basis replace bad timesamples with random numbers drawn from gaussian

    INPUTS:

    data : (array-like) filterbank data
    nchans : (int) number of filterbank channels
    sig: (float) standard deviations away from mean to clip after

    RETURNS:

    data : (array-like) rfi-clipped data
    """

    #data
    data = data
    nchans = nchans

    #channel-by-channel, rescale data
    #each channel will have a mean 0 and a standard deviation 1

    for j in range(nchans):
        #data
        channel = data[j]
        #get mean, standard deviation
        mean,std,mask = Median_clip(channel,sig,max_iter=5,full_output=True,xtol=0.5)#edit: max_iter 10>5
        if std==0.0:
            #print 'WARNING: Channel', j, 'mean, std ',mean,std
            channel-=mean
        else:
            channel-=mean
        #divide by std
            channel/=std

    #make timeseries
    timeseries=data.sum(axis=0)

    #get median and std of the timeseries
    # the sigma for median clipping is currently the same as the rescaling
    
    med,std,mask = Median_clip(timeseries,sig,max_iter=5,full_output=True,xtol=0.5) #edit: max_iter 10>5
    print 'timeseries mean, std ', med,std

    #find where timeseries data lies outside of boundaries
    #boundaries are currently 3 standard deviaions away from
    #the timeseries median
    
    minval = med-(sig*std)
    maxval = med+(sig*std)
    print 'minval, maxval ',minval,maxval
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

        #replace bad data with median of timeseries/256 and std of channel std
        channel = (toclip*np.random.normal(loc=med/np.float(nchans),scale=chan_std,size=channel.shape))+(~toclip*channel)
        #overwrite
        data[i,:]=channel    

    return data

#################################################################################

def ClipFil(in_fil,outname,outloc,bitswap,rficlip=True,clipsig=3.,toload_samps=40000):
    """
    RFI clips an input filterbank file. Outputs second, clipped file with chosen name.

    INPUTS:

    in_fil       : (str) input filterbank to clip (must be directory location and file name)
    outname      : (str) output name for clipped filterbank
    outloc       : (str) output folder for clipped filterbank
    bitswap      : (True/False boolean) if True, 8-bit input will be written out as 32-bit
                 if False, 8-bit input will be written out as 8-bit
                 and vice-versa
    rficlip      : (True/False boolean) if True, rfi sigma clipping will be applied via
                 timeseries data
    clipsig      : (float) the clipping sigma, if rficlip==True
    toload_samps : (int) the number of timesamples (spectra) to load at once while clipping.



    """
    
    print 'Input file is: ',in_fil
    fil_names=np.array([in_fil]).flatten()
    print fil_names.shape

    print 'loading filterbank\n'
    fils=[]
    for fil in fil_names:
        fils.append(fr(fil)) #store pointers to filterbank file

    print 'Calculating start mjd and samples to read and skip\n' #note, nskips should be [0] and 
                                                                 #outsamps should be the number of timesamples in the filterbank file
                                                                 #as there is only one input.

    outsamps,nskips,startTime,nchans = CombineFilUtils_FBoverlap(fils)
    nskips=np.zeros_like(nskips)
    print '    ...samples to skip (should be [0]): {0}'.format(nskips)
    print '    ...samples to read (should be fil length): {0}'.format(outsamps)

    print 'calculating data chunking information\n'
    blocksize=toload_samps
    nchunks,remainder = CombineFilUtils_FBchunking(outsamps,blocksize)
    
    print 'initialising output filterbank\n'
    fh_out,bitrate = CombineFilUtils_InitialiseOut(fil_names,
                                                   outloc,
                                                   outname,
                                                   startTime,
                                                   bitswap=bitswap)
    
    print 'output files will be {0}-bit\n'.format(bitrate)
    ##SET WHICH DATA TYPE TO WRITE AS##
    if bitrate==8:
        outdtype = np.uint8
    elif bitrate==32:
        outdtype = np.float32
    
    print 'beginning clipping\n'
    

    nfils=len(fils) #number of filterbanks to combine
    for c in range(nchunks): #loop over chunks to write out
        data=np.zeros((nchans,blocksize,nfils)) #3D array to hold data
        chunk = 0 #initialise filterbank chunk
        for i in range(nfils): #loop over telescopes (should always be 1. This is vestigal from incoherent beam code and should be cleaned up.)
            skip=int(round(nskips[i])) #number of blocks to skip
            blockstart=int(skip+(c*blocksize)) #start sample to read
            print 'Reading/Writing chunk {0}/{1}'.format(c,nchunks)
            chunk=fils[i].readBlock(blockstart,blocksize) #read chunk
            if rficlip==True: #if rfi clipping mode is on:
                print 'RFI clipping...'
                chunk=RFIclip(chunk,nchans,clipsig) #clip rfi in data chunk
            data[:,:,i]=chunk #append telescope to data
        data=Beam(data) #create beam
        if bitrate==8: #if necessary...
            data=DownSampleBits(data) #...downsample to 8-bit
        data=data.T.flatten().astype(dtype=outdtype) #reshape the data to filterbank output (low freq to high freq t1, low freq to high freq t2, ....) and recast to 32 bit float
        sppu.File.cwrite(fh_out[0], data) #write block to filterbank file
    return
