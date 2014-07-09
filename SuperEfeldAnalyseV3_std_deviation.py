#! /usr/bin/env python3.2
print("Go")
from numpy import *
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import math
from pylab import *
from scipy import *
from scipy.stats import *
import scipy as sy
from os import *
#from scipy import optimize
from scipy.optimize import curve_fit
#from numpy import polyfit
import numpy as N
import matplotlib.pyplot as plt
import matplotlib.figure
import codecs

print('*** Time-correlation/Current calculator V16 turbo ***')
print('Includes Testparticles!')
print('')


#Data
run_group_name =  'ERCONX03PaperChangeE'
Number_of_runs = 100
timesteps = 30001
temperature_begin = 0.6
temperature_end = 0.6
temperature_inc = 0.1
efield = 0.05
efield_end = 0.25
efield_inc = 0.05



#Name for output-file:
beliebig = 'NewErrors'

#Parameter:
ladung = [0.0, 2.0 / 3.0, -2.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0]
TT = [34982.0,3639.0,1079.0,454849.0,232883.0,134771.0]#    0.1 0.2 0.3 0.4  0.5 0.6 
size_of_timestep = 0.01
#Correlator_maxtime = 300
average_free_cutoff = 10.
kurzintervall=int(30001)
Anzahlxlabels = 10
Analysestart = 8000


#Arrays and Variables
Analysezeit = range(Analysestart+1,timesteps-1)
Anzahl_Analysezeitschritte = timesteps - Analysestart - 1

Positionxlabels = arange(0,kurzintervall,kurzintervall/Anzahlxlabels)
Arrayxlabels = arange(0,kurzintervall*size_of_timestep,kurzintervall*size_of_timestep/Anzahlxlabels)
kurzezeit=N.linspace(0,kurzintervall,kurzintervall)
kurzezeit_echtezeit=N.arange(0,kurzintervall*size_of_timestep,size_of_timestep)


#print('FILENAME: ' + filename_basis)
print('RUNS: ' + str(Number_of_runs))
print('timesteps: ' + str(timesteps))
print('')

#Diverses
aktuelle_ladung = 0.0
GoForAll=True
filename_basis = ''
variancestring = ''
mittelwert = N.zeros(Number_of_runs)


#Folder
mkdir('ANALYSIS_'+run_group_name+'_'+beliebig)
#chdir('ANALYSIS_'+run_group_name+'_'+beliebig)
#mkdir('Time_Histogram_'+run_group_name+'_'+beliebig)
#mkdir('Currents_'+run_group_name+'_'+beliebig)
#mkdir('Correlators_'+run_group_name+'_'+beliebig)
#chdir('../')

print('*********************')
print('New simplified correlator routine')
b = temperature_begin
while b < (temperature_end + temperature_inc):
	a = efield
	while a < (efield_end + efield_inc):
		filename_basis = run_group_name + "T" + str(b) + 'E' + str(a) + 'R'
		filename_basis = filename_basis.replace('.', '_')
		print("*****************************************")
		print("file " + filename_basis)
		print("Efield = " + str(a))
		print("Temperatur = " + str(b))
		print("*****************************************")

		aktuelle_ladung = 0.0
		ladung = [0.0, 2.0 / 3.0, -2.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0]
		maxrunnumber = Number_of_runs + 1
		Zeit_max = timesteps
		zeit = range(0, Zeit_max)
		variance_average = 0.0
		variance_std = 0.0

		#Initialize the arrays
		sum_of_corrs=N.zeros(Zeit_max)
		variance_array = N.zeros(Number_of_runs)
		corrfunction_array = N.zeros((Zeit_max, Number_of_runs)) 
		corrfunction_array_manually = N.zeros((Zeit_max, Number_of_runs)) 
		mean_corr_array = N.zeros(Zeit_max)
		std_corr_array = N.zeros(Zeit_max)
		varianzen_array= N.zeros(Zeit_max)
		gewichte_array=N.zeros(Zeit_max)
				
		stromarray = N.zeros((Zeit_max, Number_of_runs))
		mean_strom_array = N.zeros(Zeit_max)
		std_strom_array = N.zeros(Zeit_max)

		#Loop through all runs
		runnumber = 1
		File_loaded = False
		No_file_count = 0

		while runnumber < maxrunnumber:
			if runnumber==50:
				print('Run: ' + str(runnumber))
			if runnumber==100:
				print('Run: ' + str(runnumber))
			if runnumber==150:
				print('Run: ' + str(runnumber))
			if runnumber==200:
				print('Run: ' + str(runnumber))					    
			#######################################################
			#Initialize the arrays for THIS run
			vollerCorrelator = N.zeros(Zeit_max)
			stromcomplete = N.zeros(Zeit_max)
			halberCorrelator = N.zeros(Zeit_max)
			richtigerCorrelator= N.zeros(Zeit_max)
			#######################################################
			#print('Run: ' + str(runnumber))
		
			#Read the files, numbered along this scheme name1.f11
			#File is in array of strings

			filename = filename_basis + str(runnumber) + '.f11'
			zeilen = []
			try:
				z = codecs.open(filename, mode='r')
				File_loaded = True
			except:
				print('File does not exist. Skip ' + filename)
				File_loaded = False
				No_file_count += 1

			if File_loaded == True:
				for line in z.readlines():
					zeilen.append(line)
				z.close()
				Aktiv = False
				#print('Start the loop...')
				#Loop through all lines for space-averaging of the current
				for line in zeilen:
					#Data in columns
					if not line.find('Particle-Type:1')==-1:
						Aktiv = True
						ch = 1
					if not line.find('Particle-Type:2')==-1:
						Aktiv = True
						ch = 2						
					if not line.find('Particle-Type:3')==-1:
						Aktiv = True
						ch = 3
					if not line.find('Particle-Type:4')==-1:
						Aktiv = True
						ch = 4
					if not line.find('Particle-Type:5')==-1:
						Aktiv = True
						ch = 5
					if not line.find('Particle-Type:6')==-1:
						Aktiv = True
						ch = 6	
									
					if Aktiv == True:
						col = line.split('\t')
						try:
							zeitcounter = int(col[0])
							#print zeitcounter
							stromcomplete[zeitcounter] += float(col[15])*ladung[ch]
							stromarray[zeitcounter,runnumber-1] = stromcomplete[zeitcounter]
							#if zeitcounter == 20000:
							#	print ch
							#	print(zeitcounter)
							#	print stromcomplete[zeitcounter-1]
								
	#print(str(runnumber) + '\t' + str(zeitcounter) + '\t' + str(stromcomplete[zeitcounter]))
						except:
							pass
				Aktiv = False		
				#mittelwert[runnumber-1]=N.mean(stromcomplete)
				#stromcomplete -= N.mean(stromcomplete)
				
				
				##Select average-free runs TRIAL
				#if fabs(mittelwert[runnumber-1]) > average_free_cutoff:
					#No_file_count += 1
					##print('0-Wert: ' + str(runnumber)+ '  ' + str(stromcomplete[0]) + '	' + str(mittelwert[runnumber-1]))					
				#else:	
					#chdir('ANALYSIS_'+run_group_name+'_'+beliebig)
					#chdir('Time_Histogram_'+run_group_name+'_'+beliebig)
					#print('*** Plot Nx-Histogramm for this run ***')
					#plt.hist(stromcomplete,100)
					#plt.savefig(beliebig +'_HISTO-Current_' + '_'  + str(runnumber) + '.png', dpi=300, figsize=(8, 6))
					#plt.clf()
					#chdir('../../')

					##Plot Currents
					##plt.plot(kurzezeit, stromcomplete[0:1000], linestyle='-', color='r', label=('N1'))
					##plt.xlabel('Time[steps]')
					##plt.ylabel('Current, MEAN: '+str(mittelwert[runnumber-1]))
					##plt.savefig('Part_'+str(typcount)+'_'+filename_basis+'_'+str(Number_of_runs)+'_runs_'+str(Zeit_max)+'_tsteps_'+beliebig+'.eps')
					##plt.savefig('AAF-Current' + '_T_'+ str(b) +'_'+ str(runnumber) + '.png', dpi=300, figsize=(8, 6))
					##plt.clf()

					##print('Calculate the Variance manually')
					#variance = sum(pow(stromcomplete[i],2.0) for i in xrange(0, Zeit_max))/Zeit_max
					##print('Calculate the complete Correlator for this run')
					##print('----------------------------------------------')
					#variance_array[runnumber-1] = variance
					##
					##CORRELATOR:
					##
					##1) Automatic correlator
					##
					#vollerCorrelator = N.correlate(stromcomplete, stromcomplete, "full")     #/Zeit_max
					#halberCorrelator = N.array_split(vollerCorrelator, 2)
					#richtigerCorrelator = N.flipud(halberCorrelator[0])
					#N.concatenate((richtigerCorrelator.T, [0]))
					#for t in arange(Zeit_max):
						#richtigerCorrelator[t] /= float(Zeit_max - t)
						#sum_of_corrs[t] += richtigerCorrelator[t]
						#corrfunction_array[t,runnumber-1] =  richtigerCorrelator[t]
					##2) Manual Correlator
					##
					##print('Calculate the Correlator manually')
					##print('---------------------------')		
					##t=0
					##while t < Correlator_maxtime:	
					##	corrfunction_array_manually[t,runnumber-1] =  sum(stromcomplete[i]*stromcomplete[i+t] for i in arange(Zeit_max-t))/(Zeit_max-t)
					##	t = t + 1
							
			#plt.plot(zeit, stromcomplete, linestyle='-', color='r', label=('N1'))
			#plt.savefig(str(runnumber) + '_runs_' + beliebig +  '.png', dpi=300, figsize=(8, 6))
			#plt.clf()	
			runnumber = runnumber + 1
		
		
		Number_of_runs = Number_of_runs - No_file_count
		
		# Calculate the mean current and all that
		mean_strom_array = N.mean(stromarray,axis=1)
		std_strom_array = N.std(stromarray,axis=1)
		
		stdfehler_array = N.std(stromarray,axis=1)/sqrt(runnumber)
		
		
		
		
		
		#print mean_strom_array
		
		# Calculate the variances
		for i in zeit:
			varianzen_array[i]=std_strom_array[i]*std_strom_array[i]
			#print std_strom_array[i]
			if varianzen_array[i]>0:
				gewichte_array[i]=1.0/varianzen_array[i]
		#Letzter Timestep ist Null
		daten_zum_fitten = mean_strom_array[Analysestart:timesteps-1]
		fehlerdaten_zum_fitten = std_strom_array[Analysestart:timesteps-1]		
		varianzen_zum_fitten = varianzen_array[Analysestart:timesteps-1]
		gewichte_zum_fitten = gewichte_array[Analysestart:timesteps-1]
		
			
				
		# CALCULATE THE MEAN AND THE DEVIATION OF THE DATASET
		#weighted average
		Final_mean = N.average(daten_zum_fitten,weights=gewichte_zum_fitten)
		#weighted sample variance
		Weighted_Sample_variance=sum(gewichte_zum_fitten[y]*pow((daten_zum_fitten[y]-Final_mean),2.0) for y in arange(Anzahl_Analysezeitschritte))/(sum(gewichte_zum_fitten))
		Weighted_Sample_Std_Deviation = sqrt(Weighted_Sample_variance)
		
		#Inner error
		innerer_fehler = sqrt(1.0/(sum(gewichte_zum_fitten)))
		#Ausserer Fehler
		external_error = sqrt(1.0/(Anzahl_Analysezeitschritte-1.0)*1.0/sum(gewichte_zum_fitten)*sum(gewichte_zum_fitten[y]*pow((daten_zum_fitten[y]-Final_mean),2.0) for y in arange(Anzahl_Analysezeitschritte)))
		#Check
		Final_std = max(innerer_fehler,external_error)
		
		print "Weighted Mean: "
		print Final_mean
		print "internal error: "
		print innerer_fehler
		print "external error:"
		print external_error
		print "Benutzter Error:"
		print Final_std
		
		# CALCULATE THE MEAN AND THE DEVIATION OF THE DATASET
		#weighted average
		Final_mean_non_averaged = N.average(daten_zum_fitten)
		Final_mean_non_averaged_std = N.std(daten_zum_fitten)
		

		
		print "Normal Mean, weigth One"
		print Final_mean_non_averaged
		print "Normal Mean std deviation"
		print Final_mean_non_averaged_std
		print "Weighted Sample Std Deviation"
		print Weighted_Sample_Std_Deviation
		print "Weighted Sample STd ERROR"
		print Weighted_Sample_Std_Deviation/sqrt(Anzahl_Analysezeitschritte)
		print "\n\n"
		
		print(str(No_file_count) + ' files did NOT exist!!!')


		#Plot the Current	
		chdir('ANALYSIS_'+run_group_name+'_'+beliebig)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.xaxis.grid
		ax.yaxis
		figure.autolayout = True
		
		#Plot the Values:
		Meantex=r'$ \left\langle N^1\right\rangle = ' +str(round(Final_mean,6))+' \pm ' + str(round(Final_std,6)) + ' fm^{-3}$'
		#tautex=r'$'+str(round(fitParams[1]*size_of_timestep,4))+' \pm ' + str(round(sigma[1]*size_of_timestep,4)) + ' fm$'
		#ctex = r'$ c = '+ str(round(fitParams[2],4)) + ' \pm ' + str(round(sigma[2],4)) + ' fm^{-6}$'

		#tex = r'$C(0)= $' + Atex + '\n ' + r'$\tau= $' + tautex + '\n' + ctex
		ax.text(Analysestart, Final_mean*0.1, Meantex, fontsize=15, va='bottom')
		
		
		
		# Plot Values for Correlation-Functions
		#Atex=r'$'+str(round(fitParams[0],4))+' \pm ' + str(round(sigma[0],4)) + ' fm^{-6}$'
		#tautex=r'$'+str(round(fitParams[1]*size_of_timestep,4))+' \pm ' + str(round(sigma[1]*size_of_timestep,4)) + ' fm$'
		#ctex = r'$ c = '+ str(round(fitParams[2],4)) + ' \pm ' + str(round(sigma[2],4)) + ' fm^{-6}$'

		#tex = r'$C(0)= $' + Atex + '\n ' + r'$\tau= $' + tautex + '\n' + ctex
		#ax.text(fitParams[1]/2.0, fitParams[0]+fitParams[2], tex, fontsize=15, va='bottom')
		
		ax.set_xlabel(r'Time $[fm]$')
		ax.tick_params(axis='x', pad=-0.7)
		
		#ax.errorbar(x, y, yerr=[yerr_lower, 2*yerr], xerr=xerr,fmt='o', ecolor='g', capthick=2)
		plt.grid(b=True, which='major')
		#plt.plot(kurzezeit, fitFunc(kurzezeit, fitParams[0], fitParams[1], fitParams[2]))#,zeit, fitFunc(t, fitParams[0] + sigma[0], fitParams[1] - sigma[1], fitParams[2] + sigma[2]),zeit, fitFunc(t, fitParams[0] - sigma[0], fitParams[1] + sigma[1], fitParams[2] - sigma[2]))
		plt.errorbar(zeit, mean_strom_array, linestyle='-', color='r', label=('Current'), yerr=stdfehler_array, ecolor='g')

		plt.ylabel(r'$N^1(t)$')
		plt.title(r'Time evolution of $N^1$')
		#plt.legend(loc=2,shadow=True)
		#plt.yscale('log')
		#plt.xscale('log')
		#plt.ylim(0.0,5)
		#plt.xlim(0,kurzintervall)
		#xticks = Npy.arange(0.6,2.57,0.18)
		#yticks = Npy.arange(0.007,0.1,0.01)
		#plt.xaxis.set_ticks( xticks )
		#plt.yaxis.set_ticks( yticks )
		#plt.yticks(yticks,yticks)
		ax.set_xticks(Positionxlabels)
		ax.set_xticklabels(Arrayxlabels,rotation=45)
		
		#plt.ylabel(r'Correlator C(t) $[fm^{-6}]$')
		#ax.bottom = 0.25
		
		#plt.savefig('Part_'+str(typcount)+'_'+filename_basis+'_'+str(Number_of_runs)+'_runs_'+str(Zeit_max)+'_tsteps_'+beliebig+'.eps')
		plt.savefig(beliebig +'CURRENT_' + '_' + filename_basis + '_' + str(Number_of_runs) + '_runs_' + str(Zeit_max) + '_tsteps_' +  '.png', dpi=300, figsize=(10, 7))
		plt.clf()	
		
		superstring=''
		for t in zeit:
			superstring += str(t) + '\t' + str(mean_strom_array[t])+ '\t' + str(stdfehler_array[t])  + '\n'
			t = t + 1		
			
		f = codecs.open(beliebig +'STROM-Datei' + '_' + filename_basis + '_' + str(Number_of_runs) + '_runs_' + str(Zeit_max) + '_tsteps_' , 'w')
		f.write(superstring)
		f.close()	
		
		resultstring  = str(Final_mean) + '\t' + str(Final_std)+'\n'
		resultstring += str(Final_mean_non_averaged)+'\t'+str(Final_mean_non_averaged_std) + '\n'
		resultstring += str(Final_mean)+'\t'+str(Weighted_Sample_Std_Deviation) + '\n'
		resultstring += str(Final_mean)+'\t'+str(Weighted_Sample_Std_Deviation/sqrt(Anzahl_Analysezeitschritte))

		
		f = codecs.open(beliebig +'RESULTS' + '_' + filename_basis + '_' + str(Number_of_runs) + '_runs_' + str(Zeit_max) + '_tsteps_' , 'w')
		f.write(resultstring)
		f.close()			
			
			
			
			
			
			
			
			
			
			
			
		##After all runs, finish the RUN-averaging and write the data out******************************************
		#print('Do the All-Particle Plots')
		
		##Write the output-file
		#chdir('ANALYSIS_'+run_group_name+'_'+beliebig)
		#chdir('Correlators_'+run_group_name+'_'+beliebig)
		#f = codecs.open(beliebig +'Corr_' + '_' + filename_basis + '_' + str(Number_of_runs) + '_runs_' + str(Zeit_max) + '_tsteps_' , 'w')
		#f.write(superstring)
		#f.close()

		#def fitFunc(zeit, a, b, c):
    			#return a*N.exp(-1.0*zeit/b) + c
		
		#daten_zum_fitten=mean_corr_array[0:kurzintervall]
		#fehlerdaten_zum_fitten=std_corr_array[0:kurzintervall]/sqrt(Number_of_runs)
		
		#fitParams, fitCovariances = curve_fit(fitFunc, kurzezeit, daten_zum_fitten)
		#sigma = [sqrt(fitCovariances[0,0]), sqrt(fitCovariances[1,1]),sqrt(fitCovariances[2,2])]
		#print fitParams
		#print fitCovariances
		
		##Plot
		#fig = plt.figure()
		#ax = fig.add_subplot(111)
		#ax.xaxis.grid
		#ax.yaxis
		#figure.autolayout = True
		
		#Atex=r'$'+str(round(fitParams[0],4))+' \pm ' + str(round(sigma[0],4)) + ' fm^{-6}$'
		#tautex=r'$'+str(round(fitParams[1]*size_of_timestep,4))+' \pm ' + str(round(sigma[1]*size_of_timestep,4)) + ' fm$'
		#ctex = r'$ c = '+ str(round(fitParams[2],4)) + ' \pm ' + str(round(sigma[2],4)) + ' fm^{-6}$'

		#tex = r'$C(0)= $' + Atex + '\n ' + r'$\tau= $' + tautex + '\n' + ctex
		#ax.text(fitParams[1]/2.0, fitParams[0]+fitParams[2], tex, fontsize=15, va='bottom')
		
		#ax.set_xlabel(r'Time $[fm]$')
		#ax.tick_params(axis='x', pad=-0.7)
		
		##ax.errorbar(x, y, yerr=[yerr_lower, 2*yerr], xerr=xerr,fmt='o', ecolor='g', capthick=2)
		#plt.grid(b=True, which='major')
		#plt.plot(kurzezeit, fitFunc(kurzezeit, fitParams[0], fitParams[1], fitParams[2]))#,zeit, fitFunc(t, fitParams[0] + sigma[0], fitParams[1] - sigma[1], fitParams[2] + sigma[2]),zeit, fitFunc(t, fitParams[0] - sigma[0], fitParams[1] + sigma[1], fitParams[2] - sigma[2]))
		#plt.errorbar(kurzezeit, daten_zum_fitten, linestyle='-', color='r', label=('Corr'), yerr=fehlerdaten_zum_fitten, ecolor='g')

		#plt.ylabel(r'$\sigma/T$')
		#plt.title(r'Time-Correlator of $N^1$')
		##plt.legend(loc=2,shadow=True)
		##plt.yscale('log')
		##plt.xscale('log')
		##plt.ylim(0.0,5)
		#plt.xlim(0,kurzintervall)
		##xticks = Npy.arange(0.6,2.57,0.18)
		##yticks = Npy.arange(0.007,0.1,0.01)
		##plt.xaxis.set_ticks( xticks )
		##plt.yaxis.set_ticks( yticks )
		##plt.yticks(yticks,yticks)
		#ax.set_xticks(Positionxlabels)
		#ax.set_xticklabels(Arrayxlabels,rotation=45)
		
		#plt.ylabel(r'Correlator C(t) $[fm^{-6}]$')
		#ax.bottom = 0.25
		
		##plt.savefig('Part_'+str(typcount)+'_'+filename_basis+'_'+str(Number_of_runs)+'_runs_'+str(Zeit_max)+'_tsteps_'+beliebig+'.eps')
		#plt.savefig(beliebig +'Corr_' + '_' + filename_basis + '_' + str(Number_of_runs) + '_runs_' + str(Zeit_max) + '_tsteps_' +  '.png', dpi=300, figsize=(8, 6))
		#plt.clf()
		
		chdir('../')
		a = a + efield_inc
	b += temperature_inc


#f = codecs.open( beliebig+'VARIANCE' + '_' + filename_basis + '_' + str(Number_of_runs) + '_runs_' + str(Zeit_max) + '_tsteps_' , 'w')
#f.write(variancestring)
#f.close()
#chdir('../../')

    ##########################################################################################################################################################################
    ##########################################################################################################################################################################
    ##########################################################################################################################################################################



print('Finish')

