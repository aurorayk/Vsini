import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from PyAstronomy import pyasl
import bisect
import glob
import pdb
from matplotlib import cm
#from scipy import exp
import scipy.stats as ss
from scipy.optimize import curve_fit
from scipy.interpolate import spline, UnivariateSpline
from scipy import asarray as ar,exp
from scipy.signal import correlate
import emcee
from lmfit import Model
from scipy.stats import chisquare
from math import log10

#need an evenly spaced wavelength values for the rotational broadening kernel to work
def interpOntoEvenGrid(wave, flux):

    #see what the average spacing is
    medSpacing = np.median(np.diff(wave))
    #Make the new wavelength list
    newWave = np.arange(wave[0], wave[-1], medSpacing)

    #interpolate the flux and variance to the same values
    newFlux = spline(wave, flux, newWave)
    
    return newWave, newFlux

#apply the rotational broadening kernel to the evenly spaced flux and wavelength values
#returns a list of fluxes that match newWave and are artificially broadened for values from
#2.0 to 65 km/s every 4 km/s (chose a spacing of 4km/s since that is similar to IGRINS resolution) 
def makeArtificialSlowRots(newWave, newFlux, coef):

    #array of velocity values for the rotational broadening routine
    velocities = np.arange(2.0, 65.0, 4.0) 
    fluxList = []
    #loop through the velocities and add the rotationally broadened flux to a list
    for i in range(len(velocities)):
        flux = pyasl.rotBroad(newWave, newFlux, coef, velocities[i])
        fluxList.append(flux)
    
    return velocities, fluxList

def normalizeFlux(flux):
   
    normFlux = flux/np.median(flux)

    return normFlux


#define a gaussian for curve fitting 
def gaus(x, a, x0, sigma,c):

    return a*exp(-(x-x0)**2/(2*sigma**2))+c

#fit a gaussian curve to x and y coordinates
def fitGaus(x, y):
    
    mean= np.sum(x*y)/np.sum(y)
    sigma = np.sqrt(np.abs(np.sum((x - mean)**2*y) / np.sum(y)))
    popt, pcov = curve_fit(gaus, x, y, p0=[0, 0, sigma, min(y)])
    gausFit = gaus(x,*popt)
    
    return gausFit
    

#measure the full width at half max of a guassian fit to the cross correlation function
#This function should be used for lower SNR spectra 
def measureGausFWHM(RV, CC):

    #cut off regions far from the maximum because outside wings can skew gaussian
    RV = np.asarray(RV)
    CC = np.asarray(CC)
    maxRVpoint = RV[np.argmax(CC)]
    rvRegion = np.where(np.logical_and(RV>(maxRVpoint-100), RV<(maxRVpoint+100)))
    ccCut = CC[rvRegion[0]]
    rvCut = RV[rvRegion[0]]

    #fit a gaussian to the peak
    gausFit = fitGaus(rvCut, ccCut)

    #find FWHM
    maxCC = max(gausFit)
    minCC = min(gausFit)
    halfMax = ((maxCC-minCC) /2 ) + minCC
    
    firstHalfGaus = gausFit[:int(np.where(gausFit == max(gausFit)) [0])]
    secondHalfGaus = gausFit[int(np.where(gausFit == max(gausFit)) [0]):]
    lowerIndex = np.argmin(np.abs(firstHalfGaus - halfMax))
    upperIndex = np.argmin(np.abs(secondHalfGaus - halfMax))
    FWHM = rvCut[upperIndex+len(firstHalfGaus)] - rvCut[lowerIndex]
    #print(FWHM)

    #plot the gaussian to check
    #pdb.set_trace()
    plt.plot(RV, CC, color = 'r')
    plt.plot(rvCut, ccCut, color='g')
    plt.plot(rvCut, gausFit, color='k')
    plt.plot([rvCut[lowerIndex], rvCut[upperIndex+len(firstHalfGaus)]], [halfMax, halfMax])
    plt.xlabel('Radial Velocity (km/s)')
    plt.ylabel('Normalized Cross Correlation Function')
    plt.show()

    return FWHM

#measure the FWHM of the cross correlation function
#This function should be used for high SNR (~100) spectra only since it does not fit a gaussian or curve to the cc function,
#and noisy peaks in the cc function can significantly change the FWHM 
def measureFWHM(RV, CC):

    #find FWHM
    maxCC = max(CC)
    firstHalf = CC[:int(np.where(CC == maxCC) [0])]
    secondHalf = CC[int(np.where(CC == maxCC) [0]):]
    minFirstHalf = min(firstHalf)
    minSecondHalf = min(secondHalf)
    minCC = (minFirstHalf+minSecondHalf)/2
    halfMax = ((maxCC-minCC) /2 ) + minCC
    lowerIndex = np.argmin(np.abs(firstHalf - halfMax))
    upperIndex = np.argmin(np.abs(secondHalf - halfMax))
    FWHM = RV[upperIndex+len(firstHalf)] - RV[lowerIndex]
    #print(FWHM)

    return FWHM
     



#pad spectrum of the template by wraping the spectrum
def padCrossCorrelation(wavelength, flux):
    
    #see what the average spacing is
    avgSpacing = np.mean(np.diff(wavelength))
    #see what the length of the wavelengths array is
    lenWave = len(wavelength)
    #figure out the lowest wavelength value to have 3 wrapped template spectra
    lowWave = wavelength[0] - (avgSpacing*lenWave)
    highWave = wavelength[-1] + (avgSpacing*(lenWave+1))
    #make a new wavelength grid
    newWave = np.arange(lowWave, highWave, avgSpacing)
    #append 3 of the flux grids together
    newFlux = np.concatenate((flux, flux, flux))

    return newWave, newFlux

#cross correlation
def crossCorrelate(wavelength, flux, waveTemplate, fluxTemplate):

    #need the flux to be evenly spaced in log wavelength units
    logWave = np.log10(wavelength)
    logWaveTemp = np.log10(waveTemplate)
    
    #want my rv grid to be spaced every 0.1 km/s
    newLogWave = np.arange(min(logWave),max(logWave), 0.1*np.log10(np.e)/(2.998*10**5))

    #interpolate the fluxes onto that grid 
    newFlux = spline(logWave, flux, newLogWave)
    newFluxTemp = spline(logWaveTemp, fluxTemplate, newLogWave)

    #pad the template but only need ~150 km/s
    paddedWaveTemp, paddedFluxTemp = padCrossCorrelation(newLogWave, newFluxTemp)
    paddedWaveTemp = paddedWaveTemp[len(newLogWave)-2000:-(len(newLogWave)-2000)]
    paddedFluxTemp = paddedFluxTemp[len(newLogWave)-2000:-(len(newLogWave)-2000)]

    #cross correlate
    cc = np.correlate(paddedFluxTemp, newFlux, mode='valid')

    #find the corresponding radial velocity grid
    rvGrid = []
    for i in range(len(cc)):
        rv = (10**(newLogWave[-1])/10**(paddedWaveTemp[len(newLogWave)+i-1]) -1 )*2.998*10**5
        rvGrid.append(rv) 

    return rvGrid, cc




    
    
    


     

    

    

    
    
