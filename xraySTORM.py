
# coding: utf-8

# In[1]:


"""
This code (xraySTORM) was written by Dylan O'Connell (dylan.oconnell@monash.edu). 
Under the supervision of Marcus J. Kitchen and Kaye S. Morgan
It contains modular functions that build up to finding the co-ordinates of 
intensity peaks of low visibility in an image. 
"""

#import modules
import time
import os
import skimage
import skimage.io as io
import glob
import itertools
from itertools import chain
import scipy.optimize as opt
from scipy import ndimage
import skimage
import numpy as np
from PIL import Image


#Function to open a .tif image at location 'path' into a greyscale numpy array 'image'
def imit(path):
    image = io.imread(path)
    return image

#Function to convert a greyscale image 'array' into a binary image 'arr' with a given threshold 'greyvalue'
def binit(array, greyvalue):
    arr = 1*(array > greyvalue)
    return arr


#Returns flattened expression of a symmetric two-D Gaussian given parametres 
def twoD_Gaussian(xdata_tuple, height, center_x, center_y, width_, offset):
    (y, x) = xdata_tuple 
    center_x = float(center_x)
    center_y = float(center_y)    
    g = height*np.exp(-(((center_x-x)/width_)**2+((center_y-y)/width_)**2)/2) + offset
    return g.ravel()

#Given a two-D array, 'array', this function will map it the function 'twoD_Gaussian' given lists of
#'intials' &  'bounds'
def gfit(array,initials, bounds):
    y = np.linspace(0,array.shape[0],array.shape[0],dtype=int)
    x = np.linspace(0,array.shape[1],array.shape[1], dtype=int)
    y, x = np.meshgrid(y, x)
    popt, pcov = opt.curve_fit(twoD_Gaussian, (y, x), array.ravel(), p0=initials, bounds=bounds)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr


#Given a location array from scipys ndimage.label package, this function will 
#segement a small region around the approximate co-ordinates of a photon cloud localisation, 
#then use 'gfit' to map to this small region.
#It will do this for every localisation within an image. 
#it will then save the paramaters for each localisation into a csv file given a savepath
def fitablob(locationarray, originalarray, numberofblobs,imagename, savepath):
    results = np.asarray([0,0,0,0,0,0,0,0,0,0,0])
    
    #Determine name of csv file to contain results
    junk, nametif = imagename.rsplit('/',1)
    name, tif = nametif.rsplit('.',1)
    
    #iterate through each localisation in an image from a localisation array
    for i in range(1,numberofblobs+1):
        #find approximate co-ordinates for localisation
        justone = (locationarray == i)*locationarray
        firstpoint = np.argmax(justone)
        y,x = np.unravel_index(firstpoint, (locationarray.shape[0],locationarray.shape[1]))
        n = 5
        if y >n:
            if x>n:
                #find approximate region for localisation
                area = originalarray[y-n:y+n, x-n:x+n]
                maxpoint = np.argmax(area)
                guessy,guessx = np.unravel_index(maxpoint, (area.shape[0],area.shape[1]))
                newpy = y + guessy - n
                newpx = x + guessx - n

                if newpy < (locationarray.shape[0] - 3):
                    if newpx < (locationarray.shape[1] - 3):
                        m = 3
                        #center a 7 by 7 region around the maximum value in approximate region
                        area = originalarray[newpy-m+1-1:newpy+m+1, newpx-m+1-1:newpx+m+1]
                        
                        #find toal amount of adc counts in a region
                        integrate = np.sum(area)
                        
                        #give intials and bounds to for fitting function
                        initials  = originalarray[newpy,newpx], 3, 3, 0.9 ,110
                        lowers = [10, 1.5, 1.5, 0.5, 90]
                        highers = [300, 4.5,4.5, 2 ,120]
                        
                        #now fit using gifit function
                        try:
                            params, stderr = gfit(area, initials, bounds=(lowers,highers))
                        #if runtime or value error occurs move onto next localisation
                        except RuntimeError:
                            continue
                        except ValueError:
                            continue
                        #State a maximum allowed uncertainty in pixels for each x center and y center value
                        allowed  = 0.5
                        #If uncertainty is larger, move onto next localisation
                        if stderr[1] > allowed:
                            continue
                        if stderr[2] > allowed:
                            continue
                        if params is None:
                            continue

                            
                        #organise paramters before saving    
                        height, center_x, center_y, width_, offset = params
                        params[2] = params[2] + newpy - m
                        params[1]= params[1] + newpx - m
                        oneresult = list(zip(params, stderr))
                        oneresult = list(chain(*oneresult))
                        oneresult.append(integrate)
                        oneresult = np.asarray(oneresult)
                        oneresult = np.around(oneresult, decimals=3)
                        results = np.vstack((results,oneresult))
    #save csv
    np.savetxt(savepath+'/'+str(name)+".csv", results, delimiter=",", 
               header="h,u(h),x,u(x),y,u(y),width,u(width),offset,u(offset),area")

    return []

#given an image, a thresholding value 'binn', a structuring element 'structure'
#and a saving directory 'savepath'. It will find and fit each photon, saving the results in csv format  
def xraySTORM(image, binn,structure, savepath):
    #open image
    im = io.imread(image)
    #Convert to binary
    origbin = binit(im, binn)
    #Apply binary erosion-dilation
    openbig = ndimage.morphology.binary_opening(origbin, structure)
    #obtain a 'location' image, where each localisation has a different grey-scale
    location, numbers = ndimage.label(openbig)
    #fit each localisation using 'fitablob'
    fitablob(location, im, numbers,image, savepath)
    return []  




