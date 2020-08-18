# File for creating a list of all .f90 files in the version of socrates I have

import os
from os import walk

socdir = '/u/home/m/mmckinne/hab_proj/Isca/src/atmos_param/socrates/src/trunk/src'
refdir = '/u/home/m/mmckinne/hab_proj/Isca/src/extra/model/socrates/path_names_old'
newdir = '/u/home/m/mmckinne/hab_proj/Isca/src/extra/model/socrates/path_names'


# Get path names for all socrates files
f = []

for (dirpath, dirnames, filenames) in walk(socdir):
        for fname in filenames:
                fpath = dirpath+'/'+fname
                f.append(fpath)

# Remove the initial part of the path name
rmstr = '/u/home/m/mmckinne/hab_proj/Isca/src/'
fprints = []
for fpath in f:
        fprints.append(fpath[len(rmstr):]+'\n')

# Remove path names that no longer exist from original path names file, and add all new /radiance_core/ names
try:
 pn = open(refdir,'r')
except:
 os.system('cp '+newdir+' '+refdir)
 pn = open(refdir,'r')
pnlines = pn.readlines()
pn.close()

pn_new = []
for line in pnlines:
	if 'atmos_param/socrates/src/trunk/src' in line and line not in fprints:
		continue
	elif 'atmos_param/socrates/src/trunk/src/aux/interp.f' in line:
		continue
	else:
		pn_new.append(line)
                
for line in fprints:
	if 'atmos_param/socrates/src/trunk/src/radiance_core' in line and line not in pn_new:
		pn_new.append(line)

pn = open(newdir,'w+')
for line in pn_new:
        pn.write(line)
pn.close()
