import sys
import numpy as np

# time frames to ouptut colvar files for
# time ranges [STIMES, ETIMES] in ps
STIMES=[0, 2000, 4000, 6000, 8000, 0, 2000]
ETIMES=[4000, 6000, 8000, 10000, 12000, 12000, 12000]

def get_basename(pwd_name):
   """
   returns base file or dir of given path
   file name returned with file extension
   given path can have a / at the end; will be ignored
   """
   import os
   if pwd_name[-1] == '/': return os.path.basename(pwd_name[:-1])
   else: return os.path.basename(pwd_name)

def remove_duplicate_times(times, colvar):
   curr=times[0]
   ndxs=[]
   for i in range(1, len(times)):
      if times[i] <= curr: ndxs.append(i)
      else: curr = times[i]

   times = np.delete(times, ndxs)
   colvar = np.delete(colvar, ndxs, axis=0)

   return times, colvar

def main():
   nranges=len(STIMES)
   lines=[[] for i in range(nranges)]

   colvar=np.loadtxt(sys.argv[1], usecols=(0,1)) 
   times=np.around(np.loadtxt(sys.argv[1], usecols=(0,)), decimals=0)
   times, colvar = remove_duplicate_times(times, colvar)

   for i in range(nranges):
      ndxs=np.where(np.logical_and(times>=STIMES[i], times<=ETIMES[i]))[0]
      fname='%s_ps%i_%i' % (get_basename(sys.argv[1]), STIMES[i], ETIMES[i])
      with open(fname, 'w') as f:
         for j in ndxs:
            f.write('%.6f %.6f\n' % (colvar[j][0], colvar[j][1]))



if __name__ == '__main__':
   main()
