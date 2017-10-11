#!/usr/bin/env python
# coding: utf-8

from hirbin import *
from hirbin.parsers import *
#from parsers import *
import os
from os import path, mkdir
from os.path import isdir
import glob
import argparse
import sys
from multiprocessing.dummy import Pool


def parseArgs(argv):
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Perform functional annotation of assembled metagenomics sequences, using TIGRFAM or PFAM.')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('-db', dest='database_dir', help='The path to the protein family database directory (required)',required=True)
  parser.add_argument('--type', dest='type', default='nucleotide', help='Sequence type. "prot" if the reference fasta files are for amino acid sequences. "nucl" if nucleotide sequences. In the latter case the nucleotide sequences will be translated. (default:  %(default)s)')
  parser.add_argument('-o', dest='output_dir', help='The path to the output directory. If not specified, a new directory is created')
  parser.add_argument('-e','--e-value', dest='evalue_cutoff', default='1e-10', help='E-value cutoff [default = %(default)s]')
  parser.add_argument('-n', dest='n',help='number of threads [default = %(default)s]',default=1, type=int)
  parser.add_argument('-f',dest='force_output_directory',action="store_true",help='Force, if the output directory should be overwritten')
  parser.add_argument('-p','--maxOverlap',dest='max_acceptable_overlap',default=0.1,help='Max percentage of acceptable sequence overlap for bins at the same contig (if overlapping more than p percent, take the best scoring HMM-profile), default=%(default)s')
  parser.add_argument('--tentacle_format', dest='tentacle_format', action='store_true',help='Write the output in Tentacle format (default: gff format)')
  arguments=parser.parse_args(argv)
  return arguments


  
def runTranslation(sample_dict,ncpus,output_dir):
  ''' Translate the nucleotide sequences in the assembly w.r.t all 6 reading frames using transeq, using n cpu:s'''
  p=Pool(ncpus)
  arglist=[path+' '+output_dir+'/protseq/'+name+'.pep' for (name,path) in sample_dict.iteritems()]
  p.map(translateOneSample,arglist)
  return(output_dir+'/protseq/')
  
def translateOneSample(argument):
  '''Subroutine for translating one sample on one cpu'''
  os.system('transeq --frame 6 ' + argument)
  
def runHMMer(sample_dict,ncpus,output_dir,database_dir,evalue_cutoff,protseq_dir,type):
  ''' Run HMMer to perform functional annotation for all samples using n cpu:s'''
  p=Pool(ncpus)
  if type.lower().startswith('prot'):
    arglist=['--domtblout ' + output_dir +'/hmmeroutput/' + name + '.hmmout -E ' + evalue_cutoff + ' ' + database_dir + ' ' +  path for (name,path) in sample_dict.iteritems()]
  else:
    arglist=['--domtblout ' + output_dir +'/hmmeroutput/' + name + '.hmmout -E ' + evalue_cutoff + ' ' + database_dir + ' ' +  protseq_dir + '/'+ name +'.pep' for name in sample_dict.keys()]
  p.map(functionalAnnotationOneSample,arglist)

def functionalAnnotationOneSample(argument):
  ''' Run HMMer for one sample using 1 cpu'''
  os.system('hmmsearch '+ argument +'> /dev/null') #dispose unwanted output to /dev/null

def runConvertCoord(sample_dict,output_directory,protseq_dir,max_acceptable_overlap,tentacle_format,type):
  '''Convert coordinates from protein to nucleotides and save in a file for annotation'''
  if type.lower().startswith('prot'):
    for (name,path) in sample_dict.iteritems():
      if tentacle_format:
        outputfilename=output_directory+'/'+name+'.tab'
      else:
        outputfilename=output_directory+'/'+name+'.gff'
      convert_coordinates(path,output_directory+'/hmmeroutput/' + name + '.hmmout',outputfilename,tentacle_format,max_acceptable_overlap)
  else:
    for name in sample_dict:
      if tentacle_format:
        outputfilename=output_directory+'/'+name+'.tab'
      else:
        outputfilename=output_directory+'/'+name+'.gff'
      convert_coordinates(protseq_dir+'/'+name+'.pep',output_directory+'/hmmeroutput/' + name + '.hmmout',outputfilename,tentacle_format,max_acceptable_overlap)


def main(mapping_file,database_dir,output_directory,type,ncpus,evalue_cutoff,force,max_acceptable_overlap,tentacle_format):
  metadata=Hirbin_run(output_directory)
  metadata.readMetadata(mapping_file)
  output_directory=metadata.createOutputDirectory(output_directory)
  metadata.output_directory=output_directory
  if type.lower().startswith("nucl"):
    #if protein sequences are missing, run translation
    try:
      mkdir(output_directory+'/protseq/')
    except OSError as e:
      if not force:
        print "Output directory already exists, you can use an already existing output directory by including the flag -f"
        raise
    protseq_dir=runTranslation(metadata.reference,ncpus,output_directory)
    protseq_dir=output_directory+'/protseq/'
  if type.lower().startswith('prot'):
    protseq_dir=""
  try:
    mkdir(output_directory+'/hmmeroutput/')
  except OSError as e:
    if not force:
      print "Output directory already exists, you can use an already existing output directory by including the flag -f"
      raise
  runHMMer(metadata.reference,ncpus,output_directory,database_dir,evalue_cutoff,protseq_dir,type)
  runConvertCoord(metadata.reference,output_directory,protseq_dir,max_acceptable_overlap,tentacle_format,type)
  
  #perform functional annotation using HMMer

if __name__=='__main__':
  arguments=parseArgs(sys.argv[1:])
  main(arguments.mapping_file,arguments.database_dir,arguments.output_dir,arguments.type,arguments.n,arguments.evalue_cutoff,arguments.force_output_directory,arguments.max_acceptable_overlap,arguments.tentacle_format)
