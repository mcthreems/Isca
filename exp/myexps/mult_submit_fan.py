# Python code to submit multiple jobs

import os, subprocess, time
from math import floor

def get_pstr(param,pid):
 if pid == 'bd':
  if param < 1:
   pout = str(int(param*100))+'cm'
  else:
   pout = str(int(param))+'m'
 elif pid == 'es0':
  if not isinstance(param, int):
   nstr = str(param)
   dloc = nstr.index('.')
   pout = nstr[:dloc]+'p'+nstr[dloc+1:]
  else:
   pout = str(int(param))
 else:
  if not isinstance(param, int):
   pfloor = floor(param)
   pout = str(pfloor)+'p'+str(int((param-pfloor)*100))+'d'
  else:
   pout = str(int(param))+'d'

 return pout

pref = 'fandry'
betas = [] #[1.e-3,3.e-3,1.e-2,3.e-2,1.e-1,3.e-1,5.e-1,7.e-1]
es0s = [1.5,1.75] #[1.e-3,3.e-3,1.e-2,3.e-2,1.e-1,3.e-1,5.e-1,7.e-1,1,1.25,1.5,1.75,2]

sublim = 4

for beta in betas:
 otxt = str(subprocess.check_output(['myjobs']))
 nsubs = otxt.count('  r  ') + otxt.count('  qw  ')
 while nsubs >= sublim:
  time.sleep(500)
  otxt = str(subprocess.check_output(['myjobs']))
  nsubs = otxt.count('  r  ') + otxt.count('  qw  ')

 fname = pref+'_beta_'+get_pstr(beta,'es0')+'.py'
 os.system('python '+fname+' sub')

for es0 in es0s:
 otxt = str(subprocess.check_output(['myjobs']))
 nsubs = otxt.count('  r  ') + otxt.count('  qw  ')
 while nsubs >= sublim:
  time.sleep(500)
  otxt = str(subprocess.check_output(['myjobs']))
  nsubs = otxt.count('  r  ') + otxt.count('  qw  ')

 fname = pref+'_es0_'+get_pstr(es0,'es0')+'.py'
 os.system('python '+fname+' sub')
    
