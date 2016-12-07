#!/usr/bin/env python
# coding: utf-8
def getClusterStruct(filename):
    '''Function for reading one cluster structure for one PFAM domain and save it into dictonary
    filename, file name of cluster file
    returns a dictonary with contigs/genes belonging to each cluster.
    '''
    domainID=filename.split('/')[-1].split('.')[0]
    f=open(filename)
    clusters={}
    sum=0
    for line in f:
        if line[0] in ["S","H"] :
           line=line.split("\t")
           clusterID=int(line[1])
           target=line[8].split()[0]
           if clusterID not in clusters:
              clusters[clusterID]=[]
              clusters[clusterID].append(target)
           else:
                clusters[clusterID].append(target)
       
    f.close()
    return clusters
    