#!/usr/bin/env python
""" MAGeCK count module
Copyright (c) 2014 Wei Li, Han Xu, Xiaole Liu lab 
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@version: $Revision$
@author:  Wei Li 
@contact: li.david.wei AT gmail.com
"""
from __future__ import print_function

import sys
import argparse
import math
import logging
import string
from mageck2.testVisualCount import *
from mageck2.mageckCountIO import *
from mageck2.mageckCountNorm import *
from mageck2.mageckCountQC import *



def mageckcount_parseargs():
  """
  Parse arguments. Only used when mageckCount.py is executed directly.
  """
  parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
  
  parser.add_argument('-l','--list-seq',required=True,help='A file containing the list of sgRNA names, their sequences and associated genes. Support file format: csv and txt.')
  parser.add_argument('--sample-label',default='',help='Sample labels, separated by comma (,). Must be equal to the number of samples provided. Default "sample1,sample2,...".')
  parser.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  parser.add_argument('--trim-5',type=int,default=0,help='Length of trimming the 5\' of the reads. Default 0')
  parser.add_argument('--sgrna-len',type=int,default=20,help='Length of the sgRNA. Default 20')
  parser.add_argument('--count-n',action='store_true',help='Count sgRNAs with Ns. By default, sgRNAs containing N will be discarded.')
  parser.add_argument('--fastq',nargs='+',help='Sample fastq files, separated by space; use comma (,) to indicate technical replicates of the same sample. For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq" indicates two samples with 2 technical replicates for each sample.')
  parser.add_argument('--fastq-2', nargs='+',
                           help='Paired sample fastq files (or fastq.gz files), the order of which should be consistent with that in fastq option.')
  
  args=parser.parse_args()
  
  
  return args

def mageckcount_checkargs(args):
  """
  Check arguments
  Parameter:
    args
        The parsearg object
  Return value:
    genedict
        The {sgrnaid:(sgrnaseq,geneid)} from library file
  """
  if hasattr(args,'fastq') and args.fastq :
    if args.list_seq == None:
      logging.error('No library file specified.')
      sys.exit(-1)
  if hasattr(args,'fastq_2') and args.fastq_2 :
    if len(args.fastq)!=len(args.fastq_2):
      logging.error('The number of paired fastqs must be equal to the number of fastq files provided.')
      sys.exit(-1)
  if args.sample_label!='':
    nlabel=args.sample_label.split(',')
    #nfq=args.fastq.split(',')
    nfq=(args.fastq)
    if len(nlabel)!=len(nfq):
      logging.error('The number of labels ('+str(nlabel)+') must be equal to the number of fastq files provided.')
      sys.exit(-1)
    args.sample_label=args.sample_label.split(',')
  # read library file
  genenames={} # store possible gene names
  if args.list_seq is not None:
    genedict=mageckcount_checklists(args) # sgid:(seq,gene)
    for (k,v) in genedict.items():
      genenames[v[1].upper()]=0
  else:
    genedict={}
  # check count table
  if args.count_table != None:
    genenames.clear()
    # check sample label
    if args.count_table.upper().endswith('CSV'):
      splitter=','
    else:
      splitter='\t'
    nl=0
    for line in open(args.count_table):
      nl+=1
      if nl==1:
        args.sample_label=line.strip().split(splitter)[2:]
      elif nl>1:
        tgene=line.strip().split(splitter)[1]
        genenames[tgene.upper()]=0
  # check QC module
  if args.day0_label != None:
    args.day0_label=args.day0_label.split(',')
    for sx in args.day0_label:
      if sx not in args.sample_label:
        logging.error('The control label specified in --day0-label option (' + sx + ') is not among the sample label list. Please double check.')
        sys.exit(-1)
    # check the GMT file
    if args.gmt_file == None:
      args.gmt_file=os.path.join(os.path.dirname(__file__),'mageckQC.gmt')
    # check the pathway used for enrichment
    try:
      with open(args.gmt_file,'r') as f:
        firstline=f.readline().strip().split()
        targetpathway=firstline[0]
        targetpathwayanno=firstline[1]
        pathwaygene=[x.upper() for x in firstline[2:]]
    except:
      logging.error('Cannot open pathway GMT file '+args.gmt_file+' for QC. Please double check your file location and format.')
      sys.exit(-1)
    # check whether there are overlaps between pathway and library
    pathwaygene_sel=[x for x in pathwaygene if x in genenames]
    if len(pathwaygene_sel)>0:
      logging.info('' + str(len(pathwaygene_sel)) + ' out of ' + str(len(pathwaygene)) + ' genes in ' +targetpathway+' are found in the library. These genes will be used for QC.')
    else:
      logging.warning('Found 0 genes in ' + args.gmt_file + ' that apper in your screening library.')
      # sys.exit(-1)
  return genedict

# two versions of rev comp, depending on different versions of python
'''
Reverse complement
'''
if sys.version_info >= (3,1):
  trans_table=str.maketrans("ACGT","TGCA")
  def mageckcount_revcomp(x):
    return x.translate(trans_table)[::-1]
else:
  trans_table=string.maketrans("ACGT","TGCA")
  def mageckcount_revcomp(x):
    return x.translate(trans_table)[::-1]


def mageckcount_mergedict(dict0,dict1):
  '''
  Merge all items in dict1 to dict0.

  Parameters
  ----------
  dict0
    A {seq:[count_0,...,count_k]} dictionary 
  dict1
    A {seq:count} dictionary 

  Return value
  ----------
  dict0
    Will be updated as {seq:[count_0,...,count_k,count_(k+1)]} 
  '''
  nsample=0
  if len(dict0)>0:
    #nsample=len(dict0[dict0.keys()[0]])
    nsample=len(list(dict0.values())[0]) 
  for (k,v) in dict0.items():
    if k in dict1:
      v+=[dict1[k]]
    else:
      v+=[0]
  for (k,v) in dict1.items():
    if k not in dict0:
      if nsample>0:
        dict0[k]=[0]*nsample
      else:
        dict0[k]=[]
      dict0[k]+=[v]
  # return dict0




def mageckcount_mergeumidict(dict0,dict1):
  '''
  Merge all items in dict1 to dict0 for umi counts.

  Parameters
  ----------
  dict0
    A  {seq:{umi:[count_0,...,count_k]}} 
  dict1
    A  {seq:{umi:count}} (dict)

  Return value
  ----------
  dict0
    Will be updated as {seq:{umi:[count_0,...,count_k,count_(k+1)]}} 
  '''
  nsample=0
  if len(dict0)>0:
      umidict=list(dict0.values())[0]
      nsample=len(list(umidict.values())[0]) 
  for (k,umidict0) in dict0.items():
    if k in dict1:
      umidict1=dict1[k]
      # merge umi counts
      # umidict0 is {umi:[count_0,...count_k]}
      # umidict1 is {umi:count}
      mageckcount_mergedict(umidict0,umidict1)
    else:
      # add "0"s 
      for (umi,countv) in umidict0.items():
        countv+=[0]
  # for the remaining guides not in dict0, fill up the records in dict0
  for (k,umidict1) in dict1.items():
    if k not in dict0:
      dict0[k]={}
      for umi in umidict1.keys():
        if nsample>0:
          dict0[k][umi]=[[0]*nsample, umidict1[umi]]
        else:
          dict0[k][umi]=[umidict1[umi]]
  # return dict0

def mageckcount_squashdict(dict0,dict1,isdict=False):
  '''
  Squash dictionaries of technical read counts

  Parameters
  ----------
  dict0
    A {seq:count_0} dictionary 
  dict1
    A {seq:count} dictionary 
  isdict
    If set, dict0 and dict1 will become {seq:{umi:count}}

  Return value
  ----------
  dict0
    Will be updated as {seq:count_0+count]} 

  '''
  if isdict==False:
    for (k,v) in dict1.items():
      if k not in dict0:
        dict0[k]=0
      dict0[k]+=v
  else:
    for (k, umidict1) in dict1.items():
      if k not in dict0:
        dict0[k]={}
      mageckcount_squashdict(dict0[k],dict1[k],isdict=False)

def mageckcount_getsamplelabel(args,datastat):
  '''
  Get a list of sample labels
  '''
  slabel=[]
  if args.fastq is not None:
    allfastq=args.fastq
    nsample=len(allfastq)
    slabel=[datastat[f.split(',')[0]]['label'] for f in allfastq]
  elif args.count_table is not None:
    nsample=len(datastat)
    slabel=['unknown']*nsample
    for (label,dsitem) in datastat.items():
      slabel[dsitem['sampleindex']]=label
    #slabel=datastat.keys()
  return slabel

def mageckcount_printdict(dict0,args,ofile,ounmappedfile,sgdict,datastat,sep='\t'):
  '''
  Write the table count to file
  '''

  slabel=mageckcount_getsamplelabel(args,datastat)
  # print header
  print('sgRNA'+sep+'Gene'+sep+sep.join(slabel),file=ofile)
  # print items
  if len(sgdict)==0:
    for (k,v) in dict0.items():
      print(sep.join([k,'None']+[str(x) for x in v]),file=ofile)
  else:
    for (k,v) in dict0.items():
      if k not in sgdict: # only print those in the genedict
        if ounmappedfile != None:
          print(sep.join([k,k]+[str(x) for x in v]),file=ounmappedfile)
        continue
      sx=sgdict[k]
      print(sep.join([sx[0],sx[1]]+[str(x) for x in v]),file=ofile)
    # print the remaining counts, fill with 0
    for (k,v) in sgdict.items():
      if k not in dict0:
        print(sep.join([v[0],v[1]]+["0"]*nsample),file=ofile)

def mageckcount_printumidict(dict0,args,ofile,ounmappedfile,sgdict,datastat,sep='\t'):
  '''
  Write the table count to file
  '''

  slabel=mageckcount_getsamplelabel(args,datastat)
  # print header
  print('sgRNA_umi'+sep+'Gene'+sep+sep.join(slabel),file=ofile)
  # print items
  if len(sgdict)==0:
    for (k,umidict) in dict0.items():
      for (umi,v) in umidict.items():
        print(sep.join([k+'_'+umi,'None']+[str(x) for x in v]),file=ofile)
  else:
    for (k,umidict) in dict0.items():
      if k not in sgdict: # only print those in the genedict
        if ounmappedfile != None:
          for (umi,v) in umidict.items():
            print(sep.join([k+'_'+umi,k]+[str(x) for x in v]),file=ounmappedfile)
        continue
      sx=sgdict[k]
      for (umi,v) in umidict.items():
        print(sep.join( [sx[0]+'_'+umi,sx[1]] + [str(x) for x in v]),file=ofile)
    # print the remaining counts, fill with 0
    for (k,v) in sgdict.items():
      if k not in dict0:
        for (umi,v) in umidict.items():
          print(sep.join([v[0],v[1]]+["0"]*nsample),file=ofile)



def mageck_printdict(dict0,args,sgdict,sampledict,sampleids):
  """Write the normalized read counts to file
  
  Parameters
  ----------
  dict0 : dict
    a {sgRNA: [read counts]} structure
  args : class
    a argparse class
  sgdict: dict
    a {sgrna:gene} dictionary
  sampledict: dict
    a {sample name: index} dict
  sampleids: list
    a list of sample index. Should include control+treatment
  
  """
  # print header
  # print items
  dfmt="{:.5g}"
  ofile=open(args.output_prefix+'.normalized.txt','w')
  # headers
  mapres_list=['']*len(sampledict)
  for (k,v) in sampledict.items():
    mapres_list[v]=k
  if len(sampledict)>0:
    cntheader=[mapres_list[x] for x in sampleids]
  else:
    cntheader=None
  logging.info('Writing normalized read counts to '+args.output_prefix+'.normalized.txt')
  if cntheader !=None:
    print('sgRNA\tGene\t'+'\t'.join(cntheader),file=ofile)
  if len(sgdict)==0:
    for (k,v) in dict0.items():
      print(k+'\t'+'None'+'\t'+'\t'.join([str(x) for x in v]),file=ofile)
  else:
    for (k,v) in dict0.items():
      if k not in sgdict: # only print those in the genedict
        logging.warning(k+' not in the sgRNA list')
        continue
      print('\t'.join([k,sgdict[k]])+'\t'+'\t'.join([str(x) for x in v]),file=ofile)
    # print the remaining counts, fill with 0
  ofile.close()




def mageckcount_checklists(args):
  """
  Read sgRNA library file
  Parameters:
    args
        The argparse object
  Return value:
    genedict
        The {sgRNAid:(seq,geneid)} object
  
  Including sgRNAs and associated sequences and lists in csv or txt file
  format: sgRNAid  seq  geneid
  """
  genedict={}
  hascsv=False
  if args.list_seq.upper().endswith('CSV'):
    hascsv=True
  n=0
  seqdict={}
  ndup=0
  for line in open(args.list_seq):
    if hascsv:
      field=line.strip().split(',')
    else:
      field=line.strip().split('\t')
    n+=1
    if field[0] in genedict:
      logging.warning('Duplicated sgRNA label '+field[0]+' in line '+str(n)+'. Skip this record.')
      continue
    if len(field)<3:
      logging.warning('Not enough field in line '+str(n)+'. Skip this record.')
      continue
    sgrnaseq=field[1].upper()
    if n==1:
      import re
      if re.search('[^ATCG]',sgrnaseq) is not None:
        logging.info('Header line of the library file detected; skip the first line ...')
        continue
    if hasattr(args,'reverse_complement') and args.reverse_complement:
      sgrnaseq=mageckcount_revcomp(sgrnaseq)
    if sgrnaseq in seqdict:
      # logging.warning('Duplicated sgRNA sequence '+field[1]+' in line '+str(n)+'. Skip this record.')
      ndup+=1
      continue
    genedict[field[0]]=(sgrnaseq,field[2])
    seqdict[sgrnaseq]=1
  logging.info('Loading '+str(len(genedict))+' predefined sgRNAs.')
  logging.warning('There are '+str(ndup)+' sgRNAs with duplicated sequences.')
  return genedict

def mageckcount_processfastq(args,genedict,sgdict):
  """
  Main entry for fastq processing
  """
  # Check whether pair-end fastq files are provided.
  paired=False
  if args.fastq_2 != None :
    paired=True
    pairedfq=[[z for z in x.split(',')] for x in args.fastq_2]
  # listfq=args.fastq.split(',')
  listfq=[[z for z in x.split(',')] for x in args.fastq]
  nsample=len(listfq)
  # create QC statistics dictionary
  datastat={}
  # check labels
  alllabel=args.sample_label
  if alllabel=='':
    slabel=['sample'+str(x) for x in range(1,nsample+1)]
  else:
    # slabel=alllabel.split(',')
    slabel=alllabel
  for i in range(nsample):
    for fi in listfq[i]:
      datastat[fi]={}
      datastat[fi]['label']=slabel[i]
  alldict={}
  alldict_umi={}

  adjust = True
  if paired and args.count_pair.upper() == 'TRUE':
    adjust = False
  i=0
  # go through the fastq files
  for filenamelist in listfq:
    dict0={}
    dict0_umi={}
    j = 0
    for filename in filenamelist: # technical replicates; should be merged together
      if paired:
        pairedfile=pairedfq[i][j]
      else:
        pairedfile=None
      dict00={}
      dict00_umi={}
      if filename.upper().endswith('BAM'):
        mageckcount_processonefile_bam(filename,args,dict00,sgdict,datastat[filename])
      elif filename.upper().endswith('SAM'):
        mageckcount_processonefile_sam(filename,args,dict00,sgdict,datastat[filename])
      else:
        mageckcount_processonefile(filename,args,dict00,dict00_umi,sgdict,datastat[filename], pairedfile, adjust=adjust)
      # squash technical counts
      mageckcount_squashdict(dict0,dict00)
      if args.umi != 'none':
        mageckcount_squashdict(dict0_umi,dict00_umi)
      j=j+1
    i=i+1
    mageckcount_mergedict(alldict,dict0)
    if args.umi != 'none':
      mageckcount_mergeumidict(alldict_umi,dict0_umi)
  
  # write to file
  ofilel=open(args.output_prefix+'.count.txt','w')
  if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
    ounmappedfilel=open(args.output_prefix+'.unmapped.txt','w')
  else:
    ounmappedfilel=None

  mageckcount_printdict(alldict,args,ofilel,ounmappedfilel,sgdict,datastat)
  
  ofilel.close()
  if hasattr(args,'unmapped_to_file') and args.unmapped_to_file:
    ounmappedfilel.close()

  # write umi counts
  if args.umi != 'none':
    ofilel=open(args.output_prefix+'.umi_count.txt','w')
    mageckcount_printumidict(alldict_umi,args,ofilel,None,sgdict,datastat)
    ofilel.close()
  # write the median normalized read counts to csv file
  allmappeddict=None
  allmappeddict_umi=None
  if len(sgdict)>0:
    allmappeddict={k:v for (k,v) in alldict.items() if k in sgdict} # only keep those with known sgRNAs
    allmappeddict_umi={k:v for (k,v) in alldict_umi.items() if k in sgdict} # only keep those with known sgRNAs
  else:
    allmappeddict=alldict
    allmappeddict_umi=alldict_umi
  return (allmappeddict,datastat)





def getcounttablefromfile(filename):
  """
  read count table from file
  Returns:
  ---------------
  x: dict
    {sgrna:[read counts]} 
  y: dict
    {sgrna:gene}
  z: dict
    z={sample_id:index}
  """
  gtab={}
  mapptab={}
  sampleids={}
  nline=0
  nfield=-1
  # if it is CSV file
  hascsv=False
  if filename.upper().endswith('.CSV'):
    hascsv=True
  logging.info('Loading count table from '+filename+' ')
  for line in open(filename):
    nline+=1
    if nline % 100000 == 1:
      logging.info('Processing '+str(nline)+' lines..')
    try:
      if hascsv==False:
        field=line.strip().split('\t')
      else:
        field=line.strip().split(',')
      if len(field)<3:
        logging.warning('Line '+str(nline)+' of the read count table has fewer than 3 columns. Skip this line ...')
      sgid=field[0]
      geneid=field[1]
      if len(geneid)==0:
        logging.error('Error: empty gene ID is not allowed at line '+str(nline)+ '. Please double-check your read count table file.')
        sys.exit(-1)
      # check if duplicate sgRNA IDs are detected
      if sgid in gtab:
        logging.warning('Duplicated sgRNA IDs: '+sgid+' in line '+str(nline)+'. Skip this record.')
        continue
      sgrecs=[float(x) for x in field[2:]]
      # check the number of fields
      if nfield!=-1 and len(sgrecs)!=nfield:
        logging.error('Error: incorrect number of dimensions in line '+str(nline)+'. Please double-check your read count table file.')
        sys.exit(-1)
      if nline==2 and len(sampleids)>0 and len(sgrecs)!=len(sampleids):
        logging.error('Error: incorrect number of dimensions in line '+str(nline)+ ' ('+str(len(sgrecs))+')'+ ' compared with the header line (' + str(len(sampleids)) + '). Please double-check your read count table file.')
        sys.exit(-1)
      nfield=len(sgrecs)
      gtab[sgid]=sgrecs
      mapptab[sgid]=geneid
    except ValueError:
      if nline!=1:
        logging.warning('Parsing error in line '+str(nline)+'. Skip this line.')
      else:
        logging.debug('Parsing error in line '+str(nline)+' (usually the header line). Skip this line.')
        ids=field[2:]
        for i in range(len(ids)):
          sampleids[ids[i]]=i
      continue
  logging.info('Loaded '+str(len(gtab))+' records.')
  return (gtab,mapptab,sampleids)

def mageckcount_checkcontrolsgrna(args,sgrna2genelist):
  """
  Check if the provided control sgRNA is correct.
  """
  if args.control_gene != None:
    controlgenelist=[line.strip() for line in open(args.control_gene)]
    logging.info(str(len(controlgenelist))+' gene(s) used as negative controls.')
    controlsglist=[]
    for (sgid, geneid) in sgrna2genelist.items():
      if geneid in controlgenelist:
        controlsglist+=[sgid]
    logging.info(str(len(controlsglist))+' sgRNAs used as negative controls.')
    if len(controlsglist)<2:
      logging.error('Not enough control sgRNAs found in the count table. Please check your control sgRNA list.')
      sys.exit(-1)
    output_file=args.output_prefix+'.control_sgrna.txt'
    output_file_hd=open(output_file,'w')
    for sgid in controlsglist:
      print(sgid,file=output_file_hd)
    output_file_hd.close()
    args.control_sgrna=output_file
  # end if
  elif args.control_sgrna != None:
    nhit=0
    controlsglist=[line.strip() for line in open(args.control_sgrna)]
    controlsglist=list(set(controlsglist)) # get the unique sgrna list 
    gene_sg_ctrl={}
    for sgid in controlsglist:
      if sgid in sgrna2genelist:
        nhit+=1
        target_gene=sgrna2genelist[sgid]
        if target_gene not in gene_sg_ctrl:
          gene_sg_ctrl[target_gene]=0
        gene_sg_ctrl[target_gene]+=1
    # calculate the number of sgRNAs each gene has in count table
    gene_sgid={}
    for (sgid, geneid) in sgrna2genelist.items():
      if geneid not in gene_sgid:
        gene_sgid[geneid]=0
      gene_sgid[geneid]+=1
    logging.info(str(nhit)+' out of '+str(len(controlsglist))+' control sgRNAs are found in count table.')
    if nhit<2:
      logging.error('Not enough control sgRNAs found in the count table. Please check your control sgRNA list.')
      sys.exit(-1)
    # check if some genes contail both negative and non-negative controls
    for gene_check in gene_sg_ctrl:
      if gene_sg_ctrl[gene_check]<gene_sgid[gene_check]:
        logging.error('Gene '+gene_check+' consists of both negative controls (' +str(gene_sg_ctrl[gene_check])+ ') and non-negative controls (' + str(gene_sgid[gene_check]-gene_sg_ctrl[gene_check])+ '). This is not allowed -- please check your negative control sgRNA list.')
        sys.exit(-1)
      

def mageckcount_processcounttable(args,genedict,sgdict):
  """
  Main entry for count table processing
  Return value:
    gtab
        A {sgRNAID:[counts]} dict
    datastat
        Data statistics
    mapptab
        A {sgRNAID:geneID} dict
  """
  (gtab,mapptab,sampleids)=getcounttablefromfile(args.count_table)
  # check the consistency between provided library and count table
  nmiss=0
  nmissinlib=0
  if len(genedict)>0:
    for sgid in gtab.keys():
      if sgid not in genedict:
        nmiss+=1
    if nmiss>0:
      logging.warning(str(nmiss)+' sgRNA IDs in the count table are not in the library. Please double check.') 
    for sgid in genedict.keys():
      if sgid not in gtab:
        nmissinlib+=1
    if nmissinlib>0:
      logging.warning(str(nmissinlib)+' sgRNA IDs in the library are not in the count table. These sgRNAs will be counted as missing sgRNAs.') 
  # get the statistics of datasets
  datastat={}
  for (sampleid, sampleindex) in sampleids.items():
    # logging.info(sampleid+':'+str(sampleindex))
    datastat[sampleid]={}
    datastat[sampleid]['label']=sampleid
    datastat[sampleid]['sampleindex']=sampleindex
    nzerosg=0
    ntotalsgcount=0
    nrdcnt=[]
    for (sgid, cntv) in gtab.items():
      cval=cntv[sampleindex]
      if cval==0:
        nzerosg+=1
      ntotalsgcount+=cval
      nrdcnt+=[math.log(cval+1.0)]
    datastat[sampleid]['mappedreads']=ntotalsgcount
    datastat[sampleid]['reads']=ntotalsgcount
    datastat[sampleid]['zerosgrnas']=nzerosg+nmissinlib
    if len(genedict)>0:
      datastat[sampleid]['totalsgrnas']=len(genedict)
    else:
      datastat[sampleid]['totalsgrnas']=len(gtab)
    datastat[sampleid]['giniindex']=mageckcount_gini(nrdcnt)
  # end for loop for samples
  
  return (gtab,datastat,mapptab)


def mageckcount_main(args):
  """
  Main entry for mageck count module
  """
  # check arguments
  genedict=mageckcount_checkargs(args) # return: {sgrnaid:(seq,geneid)} 
  # save sgRNA ID and gene name
  sgdict={} #
  for (k,v) in genedict.items():
    sgdict[v[0]]=(k,v[1]) # {seq:(sgid,gene)
  sgrna2genelist={k:v[1] for (k,v) in genedict.items()}
  if hasattr(args,'count_table') and args.count_table != None:
    # treat it as a count table
    (allmappeddict,datastat,mapptab)=mageckcount_processcounttable(args,genedict,sgdict)
    # note that the key of allmappeddict is sgRNA ID
    # if library file is provided, we need to change sgdict to make it consistent with other situations (like fastq file)
    sgdict={k:(k,v) for (k,v) in mapptab.items()}
    sgrna2genelist=mapptab
  else:
    # check the listed files: fastq/sam/bam files provided
    (allmappeddict,datastat)=mageckcount_processfastq(args,genedict,sgdict)
    # note that the key of allmappeddict is sgRNA sequence
  
  # normalize read counts
  if hasattr(args,"norm_method"):
    normmethod=args.norm_method
  else:
    normmethod="median"
  # check control sgrna
  mageckcount_checkcontrolsgrna(args,sgrna2genelist)

  ctrlsg=args.control_sgrna

  medalldict=normalizeCounts(allmappeddict,sgdict=sgdict,method=normmethod,controlsgfile=ctrlsg)
  ofilel=open(args.output_prefix+'.count_normalized.txt','w')
  mageckcount_printdict(medalldict,args,ofilel,None,sgdict,datastat,sep='\t')
  ofilel.close()
  # perform additional QCs
  if args.day0_label!= None:
    mageckcount_getQC(args,datastat,sgdict)
  # print statistics
  mageckcount_printstat(args,datastat)
  return 0




