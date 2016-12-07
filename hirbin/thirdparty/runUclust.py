#!/usr/bin/env python
# coding: utf-8
import os,subprocess
def runUclust(path,identityCutoff,usearchpath):
  print "Clustering sequences using sequence identity cutoff " +str(identityCutoff)+'.' 
  filelist=os.listdir(path)
  outpath=path+'../clust'+str(identityCutoff)+'/'
  for f in filelist:
    if not os.path.isfile(outpath+f+'.uc'):
      #if aborted continue at previous stage
      args=[usearchpath, '-cluster_fast', path+f , '-id', str(identityCutoff), '-uc', outpath+f+'.uc']
      p=subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
if __name__=='__main__':
  import sys
  runUclust(sys.argv[1],sys.argv[2])

      