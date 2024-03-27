# Python script for generating parameter sweep experiment files

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

pref = 'tam_noinf'
basefile = pref+'_bd10m_es01_rot1d.py'

bds = [10,1,0.1,0.01] #bucket depths in meters 100,40,10,1,0.1,0.01
es0s = [1,1.25,1.5,1.75,2] #es0 values 1,1.25,1.5,1.75,2  1.e-3,3.e-3,1.e-2,3.e-2,1.e-1,3.e-1,5.e-1,7.e-1
rots = [0.5,0.75,1,4,8,16] #rotation period in days 0.5,0.75,1,4,8,16

ibase = open(basefile,'r')
btext = ibase.read()
ibase.close()

for bd in bds:
 for es0 in es0s:
  for rot in rots:
   bdstr = get_pstr(bd,'bd')
   es0str = get_pstr(es0,'es0')
   rotstr = get_pstr(rot,'rot')
   genfile = pref+'_bd'+bdstr+'_es0'+es0str+'_rot'+rotstr+'.py'
   if genfile != basefile:
    eloc1 = btext.index('bucket_depth = ')+len('bucket_depth = ')
    eloc2 = btext[eloc1:].index(' ')+eloc1
    eloc3 = btext.index('es0 = ')+len('es0 = ')
    eloc4 = btext[eloc3:].index(' ')+eloc3
    eloc5 = btext.index('omega_scale = ')+len('omega_scale = ')
    eloc6 = btext[eloc5:].index(' ')+eloc5
    gtext = btext[:eloc1]+str(bd)+btext[eloc2:eloc3]+str(es0)+btext[eloc4:eloc5]+str(rot)+btext[eloc6:]
    
    igen = open(genfile,'w+')
    igen.write(gtext)
    igen.close()

   else:
    continue 

