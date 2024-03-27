# Python code to submit multiple jobs

import os, subprocess, time

pref = 'tam_fbm'
bds = ['100m','40m','10m','1m','10cm'] #['100m','40m','10m','1m']
es0s = ['1'] #['1','1p1','1p2','1p3','1p4','1p5']
rots = ['1d','4d','8d','16d'] #['1d','4d','8d','16d']

sublim = 4

for bd in bds:
 for es0 in es0s:
  for rot in rots:
   otxt = str(subprocess.check_output(['myjobs']))
   nsubs = otxt.count('  r  ') + otxt.count('  qw  ')
   while nsubs >= sublim:
    time.sleep(500)
    otxt = str(subprocess.check_output(['myjobs']))
    nsubs = otxt.count('  r  ') + otxt.count('  qw  ')

   fname = pref+'_bd'+bd+'_es0'+es0+'_rot'+rot+'.py'
   os.system('python '+fname+' sub')
    
