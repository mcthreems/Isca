# Python script for generating parameter sweep experiment files

from math import floor

def get_pstr(param,pid):
 if pid == 'bd':
  if param < 1:
   pout = str(int(param*100))+'cm'
  else:
   pout = str(int(param))+'m'
 elif pid == 'es0':
  if param - floor(param) > 0:
   pout = str(floor(param))+'p'+str(round(10*(param - floor(param))))
  else:
   pout = str(int(param))
 else:
  pout = str(int(param))+'d'

 return pout

pref = 'ras'
basefile = pref+'_bd1m_es01_rot1d.py'

bds = [100,40,10,1,0.1,0.01] #bucket depths in meters
es0s = [1,1.1,1.2,1.3,1.4,1.5] #es0 values
rots = [1,4,8,16] #rotation period in days

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

