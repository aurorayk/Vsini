# Vsini
This code can be used to measure rotational broadening (vsini) in stellar spectra. 

How the code works: To measure the rotational broadening I first cross correlate a known slowly rotating spectrum with an artificially broadened version of itself. I then cross correlate the target star and the same unbroadened slow rotator. To estimate the rotational broadening I compare the full width at half maximums of the cross correlation functions. For more details see Kesseli et al. (2018): http://adsabs.harvard.edu/abs/2018AJ....155..225K and please cite this paper if you use the code!

To use the program, you will need a spectrum of a star you would like to calculate the rotational broadening of and the spectrum of a slowly rotating star of a similar spectral type. We find that bins of 2 - 3 spectral subtypes work well. 

Instructions: 

The file labeled measureVsiniScript.py is the one you will need to edit and the file measureVsini.py just has functions that are called. Everything that needs to be edited in measureVsiniScript.py is at the very top within a comment section that says to edit things there. The things that should be edited are: 
1. path: where the spectra files are stored
2. textfile name: edit the name of a text file that has the filenames of all the spectra you want to determine vsini values for. The way the code is written this file should only have spectra that you plan to use the same slowly rotating template and limb darkening coefficient for. The reason I did it this way is because the slowest part of the code is cross correlating the artificially broadened template with itself multiple times and if all the spectra in the list use the same template and limb darkening coefficient I only need to do this part once. This means that adding more spectra to the input file won't significantly increase your runtime but having lots of input files separately will increase the time. 
3. writefile: filename of the file that is written at the end with a summary of the results of the code. 
4. limbDarkCoef: linear limb darkening coefficient values from here (http://adsabs.harvard.edu/abs/2012yCat..35460014C) 
5. templateFilename: The slowly rotating template filename


To run the code type: python measureVsiniScript.py. The runtime is about 15 mintues for 1 spectrum but adding additional spectra is about 1 minute per additional spectrum in the input list (so like an hour for 45 spectra). 


Additional notes: 
1. For high SNR spectra (mine were ~100) I found I got lower uncertainties if I did not fit a gaussian to the cross correlation function, but a few of my very fast rotators (> 35 km/s) or for any spectra that have a lot lower SNRs you want to fit a gaussian to the cross correlation peak. To do this you can just change lines 109 and 134 in the script from "measureFWHM" to "measureGausFWHM". 

2. The script is configured now to be used for IGRINS spectra but can be used for any mid to high resolution spectra


