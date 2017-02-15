#!/usr/bin/env python

import os
import sys
import gfal2
import time

def download(input , output) :

	# Instantiate gfal2
	ctx = gfal2.creat_context()

	# Set transfer parameters
	params = ctx.transfer_parameters()
	params.overwrite = True
	params.timeout = 300

	dlCounter = 0
	isOK = False

	print 'Try to Download ' + input
	while not isOK and dlCounter < 50 :
		source = input
		destination = output
		try :
			r = ctx.filecopy(params, source, destination)
			isOK = True
		except Exception, e :
			print "Download failed : %s" % str(e)
			isOK = False
			time.sleep(20)

		dlCounter = dlCounter + 1
		
	if isOK :
		print 'Download succeeded !'
	else :
		print 'Download failed !'
		sys.exit(1)


def upload(input , output) :

	# Instantiate gfal2
	ctx = gfal2.creat_context()

	# Set transfer parameters
	params = ctx.transfer_parameters()
	params.overwrite = True
	params.timeout = 300

	dlCounter = 0
	isOK = False

	print 'Try to Upload ' + input
	while not isOK and dlCounter < 50 :
		source = input
		destination = output
		try:
			r = ctx.filecopy(params, source, destination)
			isOK = True
		except Exception, e:
			print "Upload failed : %s" % str(e)
			isOK = False
			time.sleep(20)

		dlCounter = dlCounter + 1

	if isOK :
		print 'Upload succeeded !'
	else :
		print 'Upload failed !'
		sys.exit(1)



def launch(qbar, delta, d) :

	inputFileName = 'map_' + qbar + '_' + delta + '_' + d + '.root'
	inputDir = 'srm://lyogrid06.in2p3.fr/dpm/in2p3.fr/home/calice/garillot/PolyaStudies/MulResults'
	inputFilePath = inputDir + '/' + inputFileName

	download(inputFilePath , 'file:' + inputFileName)

	os.system('./FitSim ' + inputFileName)

	outputFileName = 'Fit_' + qbar + '_' + delta + '_' + d + '.root'
	outputFile = 'srm://lyogrid06.in2p3.fr/dpm/in2p3.fr/home/calice/garillot/PolyaStudies/Fits/' + outputFileName
  
	upload('file:' + 'Fit.root' , outputFile)
  
 	os.system('rm ' + 'Fit.root')
	os.system('rm ' + inputFileName)



