#### run it: python PythonSubmit.py <no of divideJobs> <textfile containing list of files> <proper eos destination with xrootd convention>

import os, sys, re, time
from ROOT import *



def main():
    myfile=open('condor_SmallNtupleMASTER.jdl')
    jdl_lines = myfile.readlines()
    myfile.close()

    myfile=open('condor_SmallNtupleMASTER.sh')
    sh_lines = myfile.readlines()
    myfile.close()
    
    myfile=open('skimNtuples_TnPMASTER.py')
    Py_lines = myfile.readlines()
    myfile.close()

    textfileName = sys.argv[2]
    mytextfile = open(textfileName)
    filelist = mytextfile.readlines()
    mytextfile.close()

    ifile = 0
      
    fileNumber=0
    for afile in filelist:
      if '/eos/uscms/' in afile:
        afile = afile.replace('/eos/uscms/','root://cmsxrootd.fnal.gov//')
      if '/eos/cms/' in afile:
        afile = afile.replace('/eos/cms/','root://cmsxrootd.fnal.gov//')
      fileNumber=fileNumber+1
      Jobs = sys.argv[1]
      divideJobs = int(Jobs)
      myChain = TChain('ggNtuplizer/EventTree')
      inputFile=afile[:-1]
      myChain.Add(inputFile)
      EntryNumber = myChain.GetEntries()
      print 'Entry number for file ', fileNumber, ' : ', EntryNumber
      NoEntry = EntryNumber 
      for njob in xrange(0, divideJobs):
	floatchunk = NoEntry/divideJobs
	chunk = int(floatchunk)
	if(njob < (divideJobs-1)):
	  start = njob*chunk
	  stop = (njob+1)*chunk
	elif(njob == (divideJobs -1)):
	  start = njob*chunk
	  stop = NoEntry
	
	myfile = open('skimNtuples_TnP'+str(ifile)+'.py','w')
        for line in Py_lines:
            if 'ZZZZ' in line:
                line = line.replace('ZZZZ', inputFile)
            if 'OOOO' in line:
                line = line.replace('OOOO', 'Output_Skim_'+str(ifile)+'.root')
            if 'XXXX' in line:
                line = line.replace('XXXX', str(start))
            if 'YYYY' in line:
                line = line.replace('YYYY', str(stop))
            myfile.write(line)
        myfile.close()
	
        myfile = open('condor_SmallNtuple'+str(ifile)+'.sh','w')
        for line in sh_lines:
            if 'QQQ' in line:
                line = line.replace('QQQ',str(ifile))
            if 'OOOO' in line:
                line = line.replace('OOOO', 'Output_Skim_'+str(ifile)+'.root')
            if 'FINALDESTINATION' in line:
	        line = line.replace('FINALDESTINATION', sys.argv[3])
            myfile.write(line)
        myfile.close()
        
        
        myfile = open('condor_SmallNtuple'+str(ifile)+'.jdl','w')
        for line in jdl_lines:
            if 'QQQ' in line:
                line = line.replace('QQQ',str(ifile))
            myfile.write(line)
        myfile.close()
        

	os.system('condor_submit '+'condor_SmallNtuple'+str(ifile)+'.jdl')
	print '<Done>', ifile
	ifile+=1

    #print "created "+str(len(filelist))+' sets of jobs. please submit them.'
    time.sleep(120)
    ### extra security before rm ###
    os.system('pwd > check.txt')
    pwdtextfile = open('check.txt')
    pwdline = pwdtextfile.readlines()
    for aline in pwdline:
      if '/uscms_data/d3' in aline:
	#os.system('rm *')
	#os.system('rm JobCondor/*')
	print 'deleted all!'
    

       

if __name__=="__main__":
    main()
