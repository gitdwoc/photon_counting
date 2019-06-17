"""
This code (imconstruct) was written by Dylan O'Connell (dylan.oconnell@monash.edu). 
Under the supervision of Marcus J. Kitchen and Kaye S. Morgan
It contains the function "imconstructor" that converts localised x-ray interaction locations csv files output from xraySTORM.py 
into numpy arrays which can then be saved as images.
In this file an example of multiprocessing being used in conjuction is shown.
"""
#import dependant Python Libraries
import numpy as np

#required arguments are:
#a csv file "csv" produced by xraySTORM.py
#an original unprocessed image "original" to get the dimensions 
#the magnification factor "axismag". each axis will be magnified by this factor, a 2 x 2 image will become 6 x 6 if mag=3
def imconstructor(csv,original,axismag):
    
    #get dimensions or original image
    y,x = original.shape
    empty = np.zeros((int((axismag*y)),int((axismag*x))))
    
    csv = np.loadtxt(csv, delimiter=',', skiprows=(0))
    xs = axismag*(csv[:,2])

    ys = axismag*(csv[:,4])
    number = range(len(xs))
    for i in number:

        try:
            empty[int(round(ys[i])),int(round(xs[i]))] = 1
        except IndexError:
            continue

    return empty


##In practice this can be very slow as it must be performed on every csv file. so it can be benificial to set up a multiprocessing loop like the example below

#import time
#from multiprocessing import Pool
#import itertools
#from itertools import chain

# mag = 4,
# result = np.zeros((2048, 2048))
##split the amount of files you will process at one time to prevent RAM from filling for 10000 images
# split = 20
# times = 501


# if __name__ == \"__main__\":
#     tic = time.time()
#     for i in range(split*times)[0::split]:

#         p = Pool()
#         r=0
## "csvlist" is the list of csvs output from xraySTORM on an image set
#         r = p.starmap(imconstructor, zip(csvlist[i:i+split], itertools.repeat(result), itertools.repeat(axismag)))
#         res = np.sum(r, axis=0)
#         result+=res
#         result[0,0] = 0

#         p.close()
#         p.join()

#     toc = time.time()
#     print('time taken: '+str(toc-tic))
#     toc = time.time()
## "Result" will be the constructed image at the end of this loop
