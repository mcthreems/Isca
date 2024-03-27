# Python code for repeating slim_dat.py endlessly to keep scratch from hitting data limit

import os, time

fname = 'slim_dat.py'

ndays = 5                   # days to run this code for
nsec = ndays*24*60*60

waitmin = 60                # minutes in between each time it runs slim_dat.py
waitsec = waitmin*60

def nokill():
 kfile = open('slim_dat_kill.txt','r')
 ktxt = kfile.read()
 kfile.close()
 if 'kill' in ktxt:
  return False
 else:
  return True

stime = time.time()

while time.time() - stime < nsec and nokill():
 time.sleep(waitsec)
 os.system('python '+fname)

print('*ding*')
kfile = open('slim_dat_kill.txt','a')
kfile.write('\ndone')
kfile.close()

