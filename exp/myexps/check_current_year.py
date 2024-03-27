# Python script for checking latest completed year in job output files

import glob
from numpy import shape

efiles = glob.glob('*.e*')

for efile in efiles:
 fid = open(efile,'r')
 ftxt = fid.read()
 fid.close()

 maxfound = False
 rnum = 1
 while not maxfound:
  if 'Run '+str(rnum)+' complete' in ftxt:
   maxrun = rnum
   rnum += 1
  else:
   maxfound = True
 
 if rnum == 1:
  maxrun = 0

 print('\n'+efile)
 print(maxrun)

print('\n')

