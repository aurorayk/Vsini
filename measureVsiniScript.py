from measureVsini import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from PyAstronomy import pyasl
import bisect
import glob
import pdb
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy import exp
import scipy.stats as ss
from matplotlib import rc

###############
# INPUT DATA THAT SHOULD BE EDITED
###############
f = open('Vsini_input_M2-M3.txt','r') #input file: list of filenames that should all use the same template and same limb darkening coefficient

path = '/Users/aurora/Desktop/BU_Research/Vsini/IGRINSdata/plp/outdata/' #path to where files are stored

writeFile = open('vsini_output_M2-M3.txt', 'w') #the name of the output file where the summary of the results will be written

limbDarkCoeff = 0.25 #linear limb darkening coefficient from Claret (2012)

templateFilename = 'SDCK_20141011_0061.spec_a0v.fits' #filename of the slowly rotating star to use as a template

################
# End of input data -- should not need to edit code after this
################

targetList = []
for line in f:
    targetList.append(line.strip())

#use 3 orders with CO bands then other orders with good SNR and exclude order with Na doublet
ordersK = [5, 6, 7, 8, 9, 11, 12, 13, 14, 15]

vsiniList = np.zeros((len(targetList), len(ordersK)))
writeFile.write('StarName  order:[20,19,18,14,13,12]  Mean  Std  3sMean  3sStd \n')

#loop through the orders first so I don't have to remake the artificially broadened ones every time
for j in range(len(ordersK)):

        #loop through and open/find vsini for each
        for i in range(len(targetList)):
            
            print(targetList[i])

            #create a list that has the vsini's for all the orders for just one fast rotator
            vsiniListTemp = []

            pathSlow = path + templateFilename[5:13] +'/'
            pathFast = path + targetList[i][5:13] + '/'

            #open the fast rotator and the slow rotator
            slowRot = fits.open(pathSlow+templateFilename)
            fastRot = fits.open(pathFast+targetList[i])

            print(ordersK[j])
            
            myOrder = ordersK[j]
            #flux and wavelength are lists of lists with each list as a separate order
            slowWave = slowRot[1].data[myOrder] #the 0 and 1 get switched depending if you use xtellcor reduced or IGRINS pipeline reduced
            slowFlux = slowRot[0].data[myOrder]
            fastWave = fastRot[1].data[myOrder]
            fastFlux = fastRot[0].data[myOrder]
            #cut off the edges of each order
            slowWave = slowWave[100:-50]
            slowFlux = slowFlux[100:-50]
            fastWave = fastWave[100:-50]
            fastFlux = fastFlux[100:-50]

            #normalize the flux and variance
            normSlowFlux = normalizeFlux(slowFlux)
            normFastFlux = normalizeFlux(fastFlux)

            #plot the spectra to check
            #plt.plot(slowWave, normSlowFlux, fastWave, normFastFlux)
            #plt.show()
            #pdb.set_trace()
            

            ################## SLOW ROTATOR ##################
            #artificially broaden the slow rotator and cross correlate it with itself
            #only do this the first time since this list has the same slow rotator
            ##################################################
            
            if i==0:

                #call the interpOntoEvenGrid funtion since to use broaden func we need even grid
                newSlowWave, newSlowFlux = interpOntoEvenGrid(slowWave, slowFlux)
                #call the makeArtificialSlowRots function
                velocities, artBroadFluxList = makeArtificialSlowRots(newSlowWave, newSlowFlux, limbDarkCoeff)

                #make lists to store the needed quantities
                artBroadFWHMList = [] #FWHM of the artificially broadened slow rotator
                rvList = [] #list of the rv values 
                ccList = [] #list of the cc values

                #loop through the artificially broadened list
                for artBroadFlux in artBroadFluxList:

                    #cross correlate
                    rv,cc = crossCorrelate(newSlowWave, artBroadFlux, slowWave, slowFlux)
                    

                    #normalize the cross correlation function since the value does not matter
                    normcc = cc / np.median(cc)

                    #plt.plot(rv, normcc, color='b', alpha=0.5)

                    #append the normalized rv and cc to the lists
                    rvList.append(rv)
                    ccList.append(normcc)

                    tempFWHM = measureFWHM(rv, normcc)
                    artBroadFWHMList.append(abs(tempFWHM))

            ################### FAST ROTATOR #################
            #cross correlate the slow rotator with the fast rotator
            #and determine the vsini of the fast rotator
            ##################################################

            #interpolate the slow rotator onto the even grid so our process is the same as before
            newSlowWave, newSlowFlux = interpOntoEvenGrid(slowWave, slowFlux)

            #cross correlate
            rvFast, ccFast = crossCorrelate(fastWave, fastFlux, slowWave, slowFlux)

            #normalize the cross correlation function 
            normCcFast = ccFast / np.median(ccFast)

            #shift the cross correlation (RV) to the same zero as the template
            maxIndex = int(np.where(ccFast == max(ccFast) ) [0])
            rvValue = rvFast[maxIndex]

            #correct for the RV
            rvShift = rvFast - rvValue

            #plot to check
            #plt.plot(rvShift, normCcFast, color='r')
            #plt.show()

            #find the full width at half max of the xcorl using measureFWHM function
            FWHM = measureFWHM(rvShift, normCcFast)
            FWHM = abs(FWHM)

            #see what velocitity the FWHM measurement gives via a 1D linear interpolation
            vsiniFWHM = np.interp(FWHM, artBroadFWHMList, velocities)
        
            #write the vsini to the list that contains all the information
            vsiniList[i][j] = vsiniFWHM

            print('Vsini: '+str(vsiniFWHM))


#write the results
for k in range(len(targetList)):
    std = np.std(vsiniList[k])
    mean = np.mean(vsiniList[k])
    outlier = np.where(vsiniList[k] > (mean+3*std))
    sigmaClippedList = []
    for h in range(len(vsiniList[k])):
        if h in outlier:
            continue
        else:
           sigmaClippedList.append(vsiniList[k][h])
           
    mean3s = np.mean(sigmaClippedList)
    std3s = np.std(sigmaClippedList)
    
    writeFile.write(targetList[k] +'  ' + str(vsiniList[k])+ '  ' + str(mean) + '  ' + str(std) + '  ' + str(mean3s) + '  ' + str(std3s) +  '\n')






