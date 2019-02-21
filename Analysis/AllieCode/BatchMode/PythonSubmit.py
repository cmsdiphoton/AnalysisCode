#### run it: python PythonSubmit.py <files per job> <starting file number> <end file number> <proper eos destination with xrootd convention> <output file number>

import os, sys, re, time
from ROOT import *

def main():
    myfile=open('condor_anaMASTER.jdl')
    jdl_lines = myfile.readlines()
    myfile.close()

    myfile=open('condor_anaMASTER.sh')
    sh_lines = myfile.readlines()
    myfile.close()
    
    myfile=open('anaMASTER.C')
    Py_lines = myfile.readlines()
    myfile.close()

    s = sys.argv[1]
    split = int(s) 

    b = sys.argv[2]
    beginNum  = int(b) 

    end = sys.argv[3]
    endNum = int(end)

    outNum = 0
    if (len(sys.argv) > 5):
        out = sys.argv[5]
        outNum = int(out)-1
    
    total = (endNum - beginNum)+1
    mod = total%split
    print "Going to process %d files, starting with file %d and going to file %d. Processing %d files per condor job. Will have %d files left to process."%(total,beginNum,endNum,split,mod)

    ifile = outNum
    for n in range(beginNum, endNum,split):
        print n
        ifile += 1

        myfile = open('ana'+str(ifile)+'.C','w')
        for line in Py_lines:
            if 'AAAA' in line:
                line = line.replace('AAAA',str(n))
            if 'ZZZZ' in line:
                end = n+split
                line = line.replace('ZZZZ',str(end))
            if 'NNNN' in line:
                line = line.replace('NNNN', str(ifile))
            myfile.write(line)
        myfile.close()

        myfile = open('condor_ana'+str(ifile)+'.sh','w')
        for line in sh_lines:
            if 'NNNN' in line:
                line = line.replace('NNNN',str(ifile))
            if 'OOOO' in line:
                line = line.replace('OOOO', 'hist_'+str(ifile)+'.root')
            if 'FINALDESTINATION' in line:
                line = line.replace('FINALDESTINATION', sys.argv[4])
            myfile.write(line)
        myfile.close()
        
        myfile = open('condor_ana'+str(ifile)+'.jdl','w')
        for line in jdl_lines:
            if 'NNNN' in line:
                line = line.replace('NNNN',str(ifile))
            myfile.write(line)
        myfile.close()
        
        os.system('condor_submit '+'condor_ana'+str(ifile)+'.jdl')
        print '<Done>', ifile

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
