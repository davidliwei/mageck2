'''
MAGeCK MLE multiprocessing
'''

from __future__ import print_function
import re
import sys
import logging
import multiprocessing
import copy
import numpy as np
import signal

from mageck2.mleem import iteratenbem

# debug
try:
  from IPython.core.debugger import Tracer
except:
  pass

def thread_p_func(dinst,args,iteratenbemargs,returndict):
  '''
  functions for multithreading
  Parameters:
    dist
        A dictionary of instances
    iteratenbemargs
        A dictionary of arguments in iteratenbem() function
        
  '''
  name = multiprocessing.current_process().name
  ngene=0
  logging.info(name+': total '+str(len(dinst))+ ' instances.')
  for (tgid,tginst) in dinst.items():
    if ngene % 1000 ==1 or args.debug:
      logging.info(name+': Calculating '+tgid+' ('+str(ngene)+') ... ')
    iteratenbem(tginst,**iteratenbemargs)
    returndict[tgid]=tginst
    ngene+=1

  
def runem_multiproc(allgenedict,args,nproc=1, argsdict={}):
  '''
  Calling iternatembem using different number of threads
  Arguments:
    allgenedict:
        a dictionary of all gene instances
    args:
        arguments
    nproc
        The number of threads
    argsdict
        Positional arguments for iteratenbem
  '''
  # separate dicts
  instdictlist=[]
  mnger=multiprocessing.Manager()
  retdict=mnger.dict()
  if nproc<=0:
    logging.error('Error: incorrect number of threads.')
    sys.exit(-1)
  else:
    ngene=0
    instdictlist=[]
    for i in range(nproc):
      instdictlist.append({})
    for (tgid,tginst) in allgenedict.items():
      nsg=tginst.nb_count.shape[1]
      nbeta1=tginst.design_mat.shape[1]-1 # the number of betas excluding the 1st beta (baseline beta)
      if nsg>=args.max_sgrnapergene_permutation:
        logging.info('Skipping gene '+tgid+' from MLE calculation since there are too many sgRNAs. To change, revise the --max-sgrnapergene-permutation option.')
        tginst.beta_estimate=np.array([0.0]*(nsg+nbeta1))
        tginst.beta_pval=np.array([1.0]*(nbeta1))
        tginst.beta_zscore=np.array([0.0]*(nbeta1))
        tginst.beta_pval_pos=np.array([1.0]*(nbeta1))
        tginst.beta_pval_neg=np.array([1.0]*(nbeta1))
        continue
      targetlistid=ngene %nproc
      instdictlist[targetlistid][tgid]=tginst
      ngene+=1
  # start jobs
  jobs=[]
  
  if True:
    for i in range(nproc):
      j=multiprocessing.Process(target=thread_p_func, name='Thread '+str(i),args=(instdictlist[i],args,argsdict,retdict))
      jobs.append(j)
      j.start()
      logging.info(j.name+' started.')
    
    for jj in jobs:
      jj.join()
      logging.info(jj.name+' completed.')
  else:
    # solve ctrl+c issue in https://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    mp_pool=multiprocessing.Pool(nproc)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
      for i in range(nproc):
        j=mp_pool.apply_async(func=thread_p_func, args=(instdictlist[i],args,argsdict,retdict))
        logging.info('Thread '+str(i)+' started.')
    except KeyboardInterrupt:
      mp_pool.terminate()
    else:
      mp_pool.close()
    mp_pool.join()
    
  
  # post processing
  logging.info('All threads completed.')
  # save the instance
  # Tracer()()
  for tgid in retdict.keys():
    tginst=retdict[tgid]
    allgenedict[tgid]=tginst


  
 

def iteratenbem_permutation(genedict,args,sg_per_gene=-1,background_genedict=None,debug=True,nround=100,allgenedict=None,removeoutliers=False,size_factor=None):
  '''
  Perform permutation test
  Parameters:
    genedict
        A dictionary of gene structures
    args
        Parameters
    sg_per_gene
        Specify the number of sgRNAs per gene during permutation. Set this number when permuting per gene group. Default is -1 (use the same as the target gene).
    background_genedict
        A background dictionary of gene structures. If none, set it the same as genedict
    debug
        Whether to print debug information
    nround
        Number of rounds for permutation
    allgenedict
        The original dict of genes (only used for permutation)
    removeoutliers
        parameters to pass to iteratenbem
    size_factor
        size factor
  '''
  if background_genedict is None:
    background_genedict=genedict
  logging.info('Start permuting '+str(nround)+' rounds ...')
  allsg=[]
  #desmat=genedict[genedict.keys()[0]].design_mat
  desmat=list(genedict.values())[0].design_mat
  nbeta1=desmat.shape[1]-1
  for (geneid, geneinst) in background_genedict.items():
    nsg=geneinst.nb_count.shape[1]
    nsample=geneinst.nb_count.shape[0]
    countmat=geneinst.nb_count.getT()
    sgitem=[(geneinst.w_estimate[i],countmat[i]) for i in range(nsg)]
    allsg+=sgitem
  logging.info('Collecting '+str(len(allsg))+' sgRNAs from '+str(len(background_genedict))+' genes.')
  #
  if allgenedict is None:
    genedictcopy=copy.deepcopy(background_genedict)
  else:
    genedictcopy=copy.deepcopy(allgenedict)
  ngene=len(genedictcopy)
  betazeros=np.zeros((nround*ngene,nbeta1))
  #
  betaz_id=0
  for nrd in range(nround):
    np.random.shuffle(allsg)
    #
    logging.info('Permuting round '+str(nrd)+' ...')
    nid=0
    # randomly assigning sgRNAs to genes
    for (geneid, geneinst) in genedictcopy.items():
      # specify the number of sgs for this gene structure
      if sg_per_gene == -1:
        nsg=geneinst.nb_count.shape[1]
      else:
        nsg=sg_per_gene
      nsample=geneinst.nb_count.shape[0]
      selitem=allsg[nid:nid+nsg]
      countmat=np.vstack([x[1] for x in selitem])
      w_es=np.array([x[0] for x in selitem])
      geneinst.nb_count=countmat.getT()
      geneinst.w_estimate=w_es
      nid+=nsg
      if nid>= len(allsg):
        nid=0
        np.random.shuffle(allsg)
    # end gene loop
    #iteratenbem(geneinst,debug=False,estimateeff=True,updateeff=False,removeoutliers=removeoutliers,size_factor=size_factor,logem=False)
    argsdict={'debug':False,'estimateeff':True,'updateeff':False,'removeoutliers':removeoutliers,'size_factor':size_factor,'logem':False}
    runem_multiproc(genedictcopy,args,nproc=args.threads,argsdict=argsdict)
    for (geneid, geneinst) in genedictcopy.items():
      nsg=geneinst.nb_count.shape[1]
      beta_es=geneinst.beta_estimate[nsg:]
      # Tracer()()
      betazeros[betaz_id,:]=beta_es
      betaz_id+=1
    # end gene loop
  # end permutation
  logging.info('Assigning p values...')
  assign_p_value_from_permuted_beta(betazeros,genedict)
  
  return betazeros

def assign_p_value_from_permuted_beta(betazeros,genedict):
  '''
  Assigning p values to each gene based on the permutated beta scores
  '''
  ncompare=betazeros.shape[0]*1.0
  for (geneid, geneinst) in genedict.items():
    nsg=geneinst.nb_count.shape[1]
    beta_es=geneinst.beta_estimate[nsg:]
    cp_u0=np.sum(betazeros>beta_es,axis=0)
    cp_u1=np.sum(betazeros<beta_es,axis=0)
    cp_ustack=np.vstack((cp_u0/ncompare,cp_u1/ncompare))
    cp_minval=np.min(cp_ustack,axis=0)
    #cp_minvec=np.array(cp_minval)[0]
    cp_minvec=cp_minval*2
    geneinst.beta_permute_pval=cp_minvec
    geneinst.beta_permute_pval_neg=cp_ustack[1]
    geneinst.beta_permute_pval_pos=cp_ustack[0]
    # Tracer()()



def iteratenbem_permutation_by_nsg(genedict,args,debug=True,size_f=None):
  '''
  Perform permutation test, grouped by the number of sgrnas per gene
  Parameters:
    genedict
        A dictionary of gene structures
    args
        Parameters
    debug
        Whether to print debug information
    size_f
        size factor
  '''
  genedict_group={}
  # set up background gene group
  bg_genedict={}
  if args.control_sgrna is not None:
    controlsglist=[line.strip() for line in open(args.control_sgrna)]
    ncontrolsg=0
    for (geneid, geneinst) in genedict.items():
      sgid=geneinst.sgrnaid
      nsginctrl=sum([1 for x in sgid if x in controlsglist])
      if nsginctrl>0:
        bg_genedict[geneid]=geneinst
        ncontrolsg+=len(sgid)
        if nsginctrl <len(sgid):
          logging.error('Gene '+geneid+' consists of both negative controls and non-negative controls. This is not allowed -- please check your negative control sgRNA list.')
          sys.exit(-1)
    if len(bg_genedict)==0:
      logging.error('Cannot find genes containing negative control sgRNA IDs.')
      sys.exit(-1)
    logging.info('Using '+str(len(bg_genedict))+' genes and '+str(ncontrolsg)+' sgRNAs as negative controls for permutation...')
  else:
    bg_genedict=genedict
   
  # assign genes to groups according to the number of sgRNAs per gene
  for (geneid, geneinst) in genedict.items():
    nsg=geneinst.nb_count.shape[1]
    if nsg not in genedict_group:
      genedict_group[nsg]={}
    genedict_group[nsg][geneid]=geneinst
  # perform permutation based on gene groups
  ngene_keys=sorted(genedict_group.keys())
  last_permuted_beta=None
  for ngene_i in range(len(ngene_keys)):
    ngene=ngene_keys[ngene_i]
    this_genedict=genedict_group[ngene]
    if args.max_sgrnapergene_permutation>=ngene:
      logging.info('Permuting groups of gene with '+str(ngene)+' sgRNAs per gene. Group progress: '+str(ngene_i+1)+'/'+str(len(genedict_group)))
      last_permuted_beta=iteratenbem_permutation(this_genedict,args,sg_per_gene=ngene,background_genedict=bg_genedict,nround=args.permutation_round,allgenedict=genedict,removeoutliers=args.remove_outliers,size_factor=size_f)
    else:
      logging.info('Groups of gene with '+str(ngene)+' sgRNAs per gene: assigning p values based on previous group results. Group progress: '+str(ngene_i+1)+'/'+str(len(genedict_group)))
      if last_permuted_beta is None:
        logging.error('No permutation data found. Please increase the value of --max-sgrnapergene-permutation.')
        sys.exit(-1)
      assign_p_value_from_permuted_beta(last_permuted_beta,this_genedict)
      
      
  
    
