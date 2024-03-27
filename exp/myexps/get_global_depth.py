# Python script for calculating the global mean water depth for a particular land strip

delphis = [0,5,15,25,35,45,55,65,75]
odepths = [60.0,48.08,36.79,27.3,18.96,12.05,6.66,2.85,0.7571] #global depths in meters 
 
def get_global_depth(inlat):
 if inlat in delphis:
  gdepth = odepths[delphis.index(inlat)]
  return gdepth
 else:
  print('Error: No Available Depth Value for Input Latitude') 
  return 0
