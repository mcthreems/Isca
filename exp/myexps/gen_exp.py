# Python script for generating parameter sweep experiment files

from math import floor

def get_pstr(param,pid):
 if pid == 'bd':
  pout = str(param*100)+'cm'
 else:
  pout = str(param)
 if '.' in pout:
  pout = pout.replace('.','p')
 if 'p0c' in pout:
  pout = pout.replace('p0c','c')
 return pout

#beta values 1.e-3,3.e-3,1.e-2,3.e-2,1.e-1,3.e-1,5.e-1,7.e-1,1
#lat values 0,5,15,25,35,45,55,65,75
#rot values 0.5,1,2,4,8,16
#es0 values 0.5,1,1.25,1.5,1.75,2,2.25,2.5
def get_pvals(bname):
 if bname == 'beta':
  bvals = [1]
 elif bname == 'lat':
  bvals = [0,65]
 elif bname == 'xi':
  bvals = [1]
 elif bname == 'rot':
  bvals = [1,2,4,8,16]
 return bvals

pref = 'iscahydro_qdef' #iscahydro_qdef iscaaqua_qdef _seas
p1 = 'lat'

if p1 == 'beta':
 refval1 = '1'
elif p1 == 'lat':
 refval1 = '5'
p2 = 'rot'
p3 = 'xi'
refval2 = '1'
refval3 = '1'
basefile = pref+'_'+p1+'_'+refval1+'_'+p2+'_'+refval2+'_'+p3+'_'+refval3+'.py'

def get_fname(pval1,pval2,pval3):
 fname = pref+'_'+p1+'_'+str(pval1)+'_'+p2+'_'+str(pval2)+'_'+p3+'_'+str(pval3)+'.py'
 return fname

#bstrs = [get_pstr(bval,bname) for bval in bvals]

pnames = [p1,p2,p3]
pvals1 = get_pvals(p1)
pvals2 = get_pvals(p2)
pvals3 = get_pvals(p3)

for i1 in range(0,len(pvals1)):
 for i2 in range(0,len(pvals2)):
  for i3 in range(0,len(pvals3)):
   genfile = get_fname(pvals1[i1],pvals2[i2],pvals3[i3])
   if genfile == basefile:
    continue

   pvals = [pvals1[i1],pvals2[i2],pvals3[i3]]

   for pname in pnames:
    ploc = pnames.index(pname)
    pval = pvals[ploc]
    
    if ploc == 0:
     ibase = open(basefile,'r')
    else:
     ibase = open(genfile,'r')
    btext = ibase.read()
    ibase.close()

    if pname == 'bd':
     ipname = 'bucket_depth'
    elif pname == 'lat':
     ipname = 'water_lat'
    elif pname == 'rot':
     ipname = 'omega_scale'
    elif pname == 'xi':
     ipname = 'es0'
    else:
     ipname = pname
    eloc1 = btext.index(ipname+' = ')+len(ipname+' = ')
    eloc2 = btext[eloc1:].index(' ')+eloc1
    gtext = btext[:eloc1]+str(pval)+btext[eloc2:]
    
    igen = open(genfile,'w+')
    igen.write(gtext)
    igen.close()

