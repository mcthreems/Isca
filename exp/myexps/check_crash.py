# Python script for checking for crashes

import glob

goal_run = 25
fnames = glob.glob('*.e*')

for fname in fnames:
 fid = open(fname,'r')
 ftext = fid.read()
 fid.close()

 if 'failed. See log for details.' in ftext:
  for ii in range(1,26):
   if 'Run '+str(ii)+' failed. See log for details.' in ftext:
    print(fname+': run '+str(ii))
    break 
 elif 'Traceback (most recent call last):' in ftext:
  print(fname+': run ?')
 elif 'Run '+str(goal_run)+' complete' not in ftext:
  print(fname+': Failed to complete to run '+str(goal_run))
 
