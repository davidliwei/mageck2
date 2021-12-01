#!/usr/bin/env python
"""MAGeCK test module
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
import math
import types
import logging

from mageck2.mageckCount import *
from mageck2.mageckMathFunc import *
from mageck2.fileOps import *
from mageck2.testVisual import *

from mageck2.fdr_calculation import *

from mageck2.cnv_normalization import *

try:
  from IPython.core.debugger import Tracer
except:
  pass

def modelmeanvar(ctable,method='edger'):
  """
  model the relation between mean and variance
  """
  # calculate the mean and variance
  tablemat=list(ctable.values())
  meanvalue=getMeans(tablemat)
  varvalue=getVars(tablemat)
  # update on version 0.5.8: calculate the distribution of meanvalue and remove guides with extremely high counts (defined as >4 std of the avg meanvalue)
  if len(meanvalue)>10:
    m_mean=getMeans([meanvalue])
    v_mean=getVars([meanvalue])
    logging.info('Detecting outliers from variance estimation.. Avg read count:'+str(m_mean)+', Var: '+str(v_mean))
    high_cutoff=m_mean[0]+4*(v_mean[0]**0.5)
  else:
    high_cutoff=max(meanvalue)*2
  ignore_points=sum([1 for x in meanvalue if x>high_cutoff])
  if ignore_points>0:
    logging.info('Skipping '+str(ignore_points)+' sgRNAs from variance calculation because of their extreme high counts (> 4* STD (mean counts) ).')
  # choose values with variance greater than mean
  meangood=[meanvalue[i] for i in range(len(meanvalue)) if meanvalue[i]<varvalue[i] and meanvalue[i]<high_cutoff ]
  vargood=[varvalue[i]-meanvalue[i] for i in range(len(varvalue)) if meanvalue[i]<varvalue[i] and meanvalue[i]<high_cutoff ]
  if len(meangood)<10:
    logging.warning('The variances between control replicates are too small. If they are technical replicates, merge these into one sample.')
    meangood=[meanvalue[i] for i in range(len(meanvalue)) ]
    vargood=[(lambda x: x if x>0.01 else 0.01 )(varvalue[i]-meanvalue[i]) for i in range(len(varvalue)) ]
  # log
  meanglog=[math.log(x+1,2) for x in meangood]
  varglog=[math.log(x+1,2) for x in vargood]
  if method=='linear':
    # least square
    (k,b)=leastsquare(meanglog,varglog,meangood)
    if k<1:
      k=1
    if b<0:
      b=0
    return (k,b)
  elif method=='edger':
    dy=varglog
    dx=[2*x for x in meanglog]
    ret=(sum(dy)-sum(dx))*1.0/len(dx)
    return ret
  else:
    return 0


def getadjustvar(coef,meanval,method='mixed'):
  """
  From the model, get the adjusted variance
  """
  if method=='linear':
    k=coef[0];b=coef[1]
    #if type(meanval) is types.FloatType:
    if isinstance(meanval,float):
      return (meanval**k)*(2**b)+meanval
    #if type(meanval) is types.ListType:
    if isinstance(meanval,list):
      return [(z**k)*(2**b)+z for z in meanval]
  elif method=='edger':
    if isinstance(meanval,float):
      return (meanval**2)*(2**coef)+meanval
    if isinstance(meanval,list):
      return [(z**2)*(2**coef)+z for z in meanval]
  elif method=='mixed':
    var1=getadjustvar(coef,meanval,method='linear')
    var2=getadjustvar(coef[2],meanval,method='edger')
    return [ (lambda x,y: x if x>y else y)(var1[i],var2[i]) for i in range(len(var1))]
  else:
    return meanval


def calculate_gene_lfc(args,lfcval,sort_id,n_lower,sgrna2genelist,destkeys,validsgrna=None,ispos=False):
  """
  Calculate gene LFC using different methods
  Parameters:
    args
        Arguments
    lfcval
        sgRNA log fold change vector
    sortid
        sgRNA sort index
    n_lower
        alpha cutoff (integer)
    sgrna2genelist
        a {sgrnaid:gene} dict
    destkeys
        a [sgrnaid] vector
    validsgrna
        a [0/1] vector to indicate whether the sgRNA is not filtered out.
    ispos
        a boolean vector to indicate whether this is a positive selection
  Return value:
    genelfc
        a {geneid:lfc} dict
  """
  if validsgrna is None:
    validsgrna=[1]*len(lfcval) 
  genesglfc={}
  ni=0
  for i in sort_id:
    ni+=1
    if validsgrna[i]==0:
      continue
    targetgene=sgrna2genelist[destkeys[i]]
    if targetgene not in genesglfc:
      genesglfc[targetgene]=[]
    if args.gene_lfc_method=='alphamean' or args.gene_lfc_method=='alphamedian':
      if ni*1.0<=n_lower:
        genesglfc[targetgene]+=[lfcval[i]]
    else:
      genesglfc[targetgene]+=[lfcval[i]]
  genelfc={}
  for (gid,vl) in genesglfc.items():
    if args.gene_lfc_method=='median' or args.gene_lfc_method=='alphamedian':
      lfc=mmedian(vl)
    elif args.gene_lfc_method=='secondbest':
      if ispos:
        vll=sorted(vl,reverse=True)
      else:
        vll=sorted(vl)
      if len(vll)>1:
        lfc=vll[1]
      else:
        lfc=0.0
    elif args.gene_lfc_method=='mean' or args.gene_lfc_method=='alphamean':
      if len(vl)>0:
        lfc=sum(vl)/len(vl)
      else:
        lfc=0.0
    else:
      lfc=0.0
    genelfc[gid]=lfc
  return genelfc

def crispr_test(tab,ctrlg,testg,varg, destfile,sgrna2genelist,args):
  """
  main function of crispr test
  Parameters:
    tab
        Read count table
    ctrlg
        Index for control samples
    testg
        Index for treatment samples
    varg 
        Index for variance estimation samples; if it's empty, use defafult variance estimation samples
    destfile
        Prefix for output file (sgrna_summary.txt)
    sgrna2genelist
        {sgrna:gene} mapping
    args
        Arguments
  Return value:
    (lowp,highp,sgrnalfc)
    lowp
        alpha cutoff for neg. selection
    highp
        alpha cutoff for pos. selection
    lower_gene_lfc
        {gene:lfc} dict. lfc is for neg. selection
    higher_gene_lfc
        {gene:lfc} dict. lfc is for pos. selection
  """
  n=len(tab)

  # control and test matrix
  tabctrl={k:[v[i] for i in range(len(v)) if i in ctrlg] for (k,v) in tab.items()}
  tabtest={k:[v[i] for i in range(len(v)) if i in testg] for (k,v) in tab.items()}

  # # write to file
  # f = open('sgrna2genelist.txt','w')
  # import csv
  # writer = csv.writer(f,delimiter='\t')
  # [writer.writerow([k,v]) for (k,v) in sgrna2genelist.items()]
  # f.close()

  # sgrnas = tabctrl.keys()
  # f = open('tabctrl.txt','w')
  # writer = csv.writer(f,delimiter='\t')
  # [writer.writerow([sgrna] + tabctrl[sgrna]) for sgrna in sgrnas]
  # f.close()
  # f = open('tabtest.txt','w')
  # writer = csv.writer(f,delimiter='\t')
  # [writer.writerow([sgrna] + tabtest[sgrna]) for sgrna in sgrnas]
  # f.close()
  #
  #
  #
  #
  # perform copy number estimation if option selected
  if args.cnv_est is not None:
    from mageck2.cnv_estimation import mageckCNVestimation
    logging.info('Performing copy number estimation...')
    mageckCNVestimation(args.cnv_est,tabctrl,tabtest,sgrna2genelist,'CNVest',args.output_prefix)
  #
  #
  #
  #
  # control matrix for mean-var estimation
  if len(varg)>1:
    tabctrlmod={k:[v[i] for i in range(len(v)) if i in varg] for (k,v) in tab.items()}
  #elif len(ctrlg)>1 and args.variance_from_all_samples==False: # more than 1 controls
  elif len(ctrlg)>1: # more than 1 controls
    tabctrlmod={k:[v[i] for i in range(len(v)) if i in ctrlg] for (k,v) in tab.items()}
  else: # only 1 control: use all the samples for estimation
    tabctrlmod={k:[v[i] for i in range(len(v)) if i in (ctrlg+testg)] for (k,v) in tab.items()}
  # treatment matrix
  tabtreatmod={k:[v[i] for i in range(len(v)) if i in (testg)] for (k,v) in tab.items()}
  # training using control samples
  model1=modelmeanvar(tabctrlmod,method='linear')
  #model2=modelmeanvar(tabctrl,method='edger')
  model=[x for x in model1];#+[model2]
  #if type(model) is types.ListType:
  if isinstance(model,list):
    logging.debug('Adjusted model: '+'\t'.join([str(x) for x in model]))
  else:
    logging.debug('Adjusted model: k='+str(model))

  tabctrl_mat=list(tabctrl.values())
  tabctrlmodel_mat=list(tabctrlmod.values())
  tabc_mean=getMeans(tabctrl_mat)
  tabcmodel_mean=getMeans(tabctrlmodel_mat)
  #
  # setup the valid sgRNA flag; note that this step has to be done before setting the zero values of tabc_mean 
  validsgrna1=[1]*n
  if hasattr(args,"remove_zero"):
    validsgrna1=[ (lambda x: 1 if x>args.remove_zero_threshold else 0)(t) for t in tabc_mean]
  # if mean of the control samples is 0: set it to greater than 0
  tabc_min=min([x for x in tabc_mean if x>0])
  tabc_mean=[ (lambda x: x if x>tabc_min else tabc_min)(t) for t in tabc_mean]
  #
  #
  #
  #
  # calculate the variance and adjusted variance
  # for consistency, tabc_var would be the final raw variance, and tabc_adjvar would be the final adjusted variance value
  if False:
    # use only control variance
    tabc_var=getVars(tabctrlmodel_mat)
    tabc_adjvar=getadjustvar(model,tabc_mean,method='linear')
  else:
    # consider both control and treatment variances
    # raw variances
    t_var_c=getVars(tabctrlmodel_mat)
    t_var_t=getVars(tabtreatmod.values())
    n_ctrl=len(ctrlg)
    n_test=len(testg)
    if len(varg)>1 or len(ctrlg)>1 : # more than 1 controls, or users specify that variance should be calculated from certain samples
      # change back to only control variances since 0.5.7
      t_var_mix=t_var_c
    else:
      # just 1 control or users specify that variance should be calculated from all samples
      if n_ctrl==1 and n_test==1:
        t_var_mix=t_var_c # older version
      else:
        frac_c=(n_ctrl-1)*1.0/(n_ctrl-1+n_test-1)
        # frac_c=1.0 
        frac_t=1.0-frac_c
        t_var_mix=[t_var_c[i]*frac_c+t_var_t[i]*frac_t for i in range(len(t_var_c))]
        logging.info('Raw variance calculation: '+str(frac_c)+' for control, '+str(frac_t)+' for treatment')
    # adjusted variances
    tabc_var=t_var_mix
    t_var_adj=getadjustvar(model,tabc_mean,method='linear')
    if False:
      # the following code is used only if you want to consider raw variances (without regression)
      # set tabc_var and tabc_adjvar
      # calculate the fraction of raw variances in terms of adjusted variance calculation
      # if n_ctrl+n_test <=2: the variance is completely depend on modeling
      # if n_ctrl+n_test >=8: the variance is completely depend on the data
      frac_raw=(n_ctrl-1+n_test-1)/6.0
      frac_raw=0.0
      if frac_raw>1.0:
        frac_raw=1.0
      if frac_raw<0.0:
        frac_raw=0.0
      logging.info('Adjusted variance calculation: '+str(frac_raw)+' for raw variance, '+str(1-frac_raw)+' for modeling')
      tabc_var=t_var_mix
      # increase the raw variance if it's smaller than the model
      tvar_mix_upper=[max(t_var_mix[i],t_var_adj[i]) for i in range(len(t_var_mix))]
      #
      # calculate the final adjusted variance, based on either a mixture of raw variance or regressed variance.
      # Currently frac_raw=0, meaning that it's just regressed variance
      tabc_adjvar=[tvar_mix_upper[i]*frac_raw+t_var_adj[i]*(1.0-frac_raw) for i in range(len(tabc_var))]
    tabc_adjvar=t_var_adj 
  #
  #
  #
  #
  # testing using tebtest
  #nt=tabtest[tabtest.keys()[0]]
  nt=list(tabtest.values())[0] 
  ttmat=list(tabtest.values())
  ttmean=getMeans(ttmat)
  # setup the valid sgRNA flag
  validsgrna2=[1]*n
  validsgrna=[1]*n
  if hasattr(args,"remove_zero"):
    validsgrna2=[ (lambda x: 1 if x>args.remove_zero_threshold else 0)(t) for t in ttmean]
    if args.remove_zero=="control":
      validsgrna=validsgrna1
    elif args.remove_zero=="treatment": 
      validsgrna=validsgrna2
    elif args.remove_zero=="any": 
      validsgrna=[validsgrna1[t]*validsgrna2[t] for t in range(n)]
    elif args.remove_zero=="both":
      validsgrna=[1-(1-validsgrna1[t])*(1-validsgrna2[t]) for t in range(n)]
    else:
      validsgrna=[1]*n
  else:
    validsgrna=[1]*n
  logging.info("Before RRA, "+str(n-sum(validsgrna))+" sgRNAs are removed with zero counts in "+args.remove_zero+" group(s).")
  if sum(validsgrna)==0:
    logging.error('No sgRNAs are left after --remove-zero filtering. Please double check with your data or --remove-zero associated parameters.')
    sys.exit(-1)
  #
  #
  #
  #
  # calculate the p value
  try:
    # for consistency, use normal p values
    tt_p_lower=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=True)
    tt_p_higher=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=False)
    #tt_p_lower=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=True)
    #tt_p_higher=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=False)
  except:
    logging.error('An error occurs while trying to compute p values. Quit..')
    sys.exit(-1)
  #
  #
  #
  #
  # calculate tt_theta, used for RRA ranking 
  if False:
    # use ttmean to calculate the pvalue
    # first, convert to standard normal distribution values
    # old method: directly calculate the theta 
    tt_theta=[(ttmean[i]-tabc_mean[i])/math.sqrt(tabc_adjvar[i]) for i in range(n)]
  if False:
    # new method 1: use logP to replace theta
    tt_p_lower_small_nz=min([x for x in tt_p_lower if x>0 ])
    tt_p_higher_small_nz=min([x for x in tt_p_lower if x>0 ])
    tt_p_lower_rpnz=[max(x,tt_p_lower_small_nz*0.1) for x in tt_p_lower]
    tt_p_higher_rpnz=[max(x,tt_p_higher_small_nz*0.1) for x in tt_p_higher]
    tt_p_lower_log=[math.log(x,10) for x in tt_p_lower_rpnz]
    tt_p_higher_log=[math.log(x,10) for x in tt_p_higher_rpnz]
    tt_theta=[ (lambda i: tt_p_lower_log[i] if ttmean[i]<tabc_mean[i] else -1*tt_p_higher_log[i])(x) for x in range(n)]
  if True:
    # new method 2: use qnorm to reversely calculate theta, only on p_lows
    logging.info('Use qnorm to reversely calculate sgRNA scores ...')
    tt_theta_0=[(ttmean[i]-tabc_mean[i])/math.sqrt(tabc_adjvar[i]) for i in range(n)]
    tt_p_lower_small_nz=min([x*0.1 for x in tt_p_lower if x*0.1>0 ])
    #print('nz:'+str(tt_p_lower_small_nz))
    tt_p_lower_rpnz=[max(x,tt_p_lower_small_nz) for x in tt_p_lower]
    tt_p_lower_log=[math.log(x) for x in tt_p_lower_rpnz] # here it's e based log
    # convert log p to theta score
    T_Q_CV=QNormConverter() 
    tt_theta_low_cvt=T_Q_CV.get_qnorm(tt_p_lower_log,islog=True)
    # replace the None value to original theta score
    tt_theta_low=[tt_theta_low_cvt[i] if tt_theta_low_cvt[i] is not None else tt_theta_0[i] for i in range(n)]
    tt_theta=[ (lambda i: tt_theta_low[i] if ttmean[i]<tabc_mean[i] else tt_theta_0[i])(x) for x in range(n)]
  #
  tt_abstheta=[math.fabs(tt_theta[i]) for i in range(n)]
  # lower_score and higher_score are used to sort sgRNAs

  tt_p_lower_score=tt_theta
  tt_p_higher_score=[-1*x for x in tt_theta]
  #
  tt_p_twosided=[ (lambda x,y: 2*x if x<y else 2*y)(tt_p_lower[i],tt_p_higher[i]) for i in range(n)]
  tt_p_fdr=pFDR(tt_p_twosided,method=args.adjust_method)
  #
  #
  #
  #
  # CNV normalization of scores
  # map sgRNA to genes
  gene_list = []
  sgrna_list = list(tabctrl.keys())
  for sgrna in sgrna_list:
    if sgrna2genelist is not None:
      gene_list.append(sgrna2genelist[sgrna])
    else:
      gene_list.append('NA')
  #
  # normalize sgRNA scores from CNV and sort according to score
  CNVnorm = False
  if (args.cnv_norm is not None and args.cell_line is not None) or args.cnv_est is not None:
    from mageck2.cnv_normalization import read_CNVdata,sgRNAscore_piecewisenorm, highestCNVgenes
    logging.info('Performing copy number normalization.')
    # reading CNV from known CNV profiles
    if args.cnv_norm is not None and args.cell_line is not None:
      (CN_arr,CN_celldict,CN_genedict) = read_CNVdata(args.cnv_norm,[args.cell_line])
      genes2correct=False
    # predicting CNV and correcting CNV
    elif args.cnv_est is not None:
      (CN_arr,CN_celldict,CN_genedict) = read_CNVdata(args.output_prefix + 'CNVestimates.txt',['CNVest'])
      genes2correct = highestCNVgenes(CN_arr,CN_genedict,percentile=98)
    
    # correcting CNV effect
    if args.cell_line in CN_celldict or 'CNVest' in CN_celldict:
      if args.cell_line in CN_celldict:
        logging.info('Normalizing by copy number with ' + args.cell_line + ' as the reference cell line.')
      elif 'CNVest' in CN_celldict:
        logging.info('Normalizing by copy number using the estimated gene copy numbers.')
      CNVnorm = True
      norm_tt_theta = sgRNAscore_piecewisenorm(tt_theta,gene_list,CN_arr,CN_genedict,selectGenes=genes2correct)
      norm_tt_abstheta=[math.fabs(norm_tt_theta[i]) for i in range(n)]
      sort_id=[i[0] for i in sorted(enumerate(norm_tt_abstheta), key=lambda x:x[1],reverse=True)]
      # replace the original values of tt_theta
      tt_theta=norm_tt_theta
      tt_abstheta=norm_tt_abstheta
    else:
      logging.warning(args.cell_line + ' is not represented in the inputted copy number variation data.')
      sort_id=[i[0] for i in sorted(enumerate(tt_abstheta), key=lambda x:x[1],reverse=True)]
  else:
    sort_id=[i[0] for i in sorted(enumerate(tt_abstheta), key=lambda x:x[1],reverse=True)]
  #
  #
  #
  #
  # calculating lower_score and higher_score to sort sgRNAs
  tt_p_lower_score=tt_theta
  tt_p_higher_score=[-1*x for x in tt_theta]
  #
  #
  # calculating sgRNA log fold change 
  sgrnalfc=[math.log(ttmean[i]+1.0,2)-math.log(tabc_mean[i]+1.0,2) for i in range(n)]
  #
  #
  #
  # write to file
  destfname=destfile+'.sgrna_summary.txt'
  destf=open(destfname,'w')
  destkeys=list(tabctrl.keys())
  dfmt="{:.5g}"
  
  #
  # output to file
  header = ['sgrna','Gene','control_count','treatment_count','control_mean','treat_mean', 'LFC', 'control_var','adj_var','score','p.low','p.high','p.twosided','FDR','high_in_treatment']
  #if CNVnorm:
  #  header += ['CNVadj_score']
  print('\t'.join(header),file=destf)
  for i in sort_id:
    # sgRNA mapping to genes?
    if sgrna2genelist is not None:
      destkeygene=sgrna2genelist[destkeys[i]]
    else:
      destkeygene='None'
    report=[destkeys[i], destkeygene, '/'.join([dfmt.format(x) for x in tabctrl_mat[i]]), '/'.join([dfmt.format(x) for x in ttmat[i]])]
    t_r=[tabc_mean[i],ttmean[i]]
    t_r+=[ sgrnalfc[i]] # log fold change
    t_r+=[tabc_var[i],tabc_adjvar[i],tt_abstheta[i],tt_p_lower[i],tt_p_higher[i],tt_p_twosided[i],tt_p_fdr[i]]
    report+=[dfmt.format(x) for x in t_r]
    report+=[ttmean[i]>tabc_mean[i]]
    #if CNVnorm:
    #  report+=[dfmt.format(norm_tt_abstheta[i])] # add CNV-adjusted sgRNA scores
    print('\t'.join([str(x) for x in report]),file=destf)
  destf.close()
  #
  #
  #
  #
  # prepare files for gene test
  if sgrna2genelist is not None:
    destfname=destfile+'.plow.txt'
    destkeys=list(tabctrl.keys())
    sort_id=[i[0] for i in sorted(enumerate(tt_p_lower_score), key=lambda x:x[1],reverse=False)]
    # output to file
    destf=open(destfname,'w')
    print('\t'.join(['sgrna','symbol','pool','p.low','prob','chosen']),file=destf)
    for i in sort_id:
      report=[destkeys[i], sgrna2genelist[destkeys[i]],'list', tt_p_lower_score[i], '1', validsgrna[i]]
      # new in 0.5.7: only print valid sgRNAs
      if validsgrna[i]==1:
        print('\t'.join([str(x) for x in report]),file=destf)
    destf.close()

    tt_p_lower_fdr=pFDR(tt_p_lower,method=args.adjust_method)
    n_lower=sum([1 for x in tt_p_lower if x <= args.gene_test_fdr_threshold])
    n_lower_valid=sum([1 for n_i in range(n) if (tt_p_lower[n_i] <= args.gene_test_fdr_threshold) and (validsgrna[n_i] == 1) ])
    #n_lower_p=n_lower*1.0/len(tt_p_lower)
    n_lower_p=n_lower_valid*1.0/sum(validsgrna) 
    logging.debug('lower test FDR cutoff: '+str(n_lower_p))
    # calculate gene lfc
    lower_gene_lfc=calculate_gene_lfc(args,sgrnalfc,sort_id,n_lower,sgrna2genelist,destkeys,validsgrna=validsgrna)
    #
    #
    #
    destfname=destfile+'.phigh.txt'
    destf=open(destfname,'w')
    destkeys=list(tabctrl.keys())
    sort_id=[i[0] for i in sorted(enumerate(tt_p_higher_score), key=lambda x:x[1],reverse=False)]
    # output to file
    print('\t'.join(['sgrna','symbol','pool','p.high','prob','chosen']),file=destf)
    for i in sort_id:
      report=[destkeys[i], sgrna2genelist[destkeys[i]],'list', tt_p_higher_score[i], '1', validsgrna[i]]
      # new in 0.5.7: only print valid sgRNAs
      if validsgrna[i]==1:
        print('\t'.join([str(x) for x in report]),file=destf)
    destf.close()
    
    tt_p_higher_fdr=pFDR(tt_p_higher,method=args.adjust_method)
    n_higher=sum([1 for x in tt_p_higher if x <= args.gene_test_fdr_threshold])
    n_higher_valid=sum([1 for n_i in range(n) if (tt_p_higher[n_i] <= args.gene_test_fdr_threshold) and (validsgrna[n_i] == 1) ])
    if n_higher>0:
      #n_higher_p=n_higher*1.0/len(tt_p_higher)
      n_higher_p=n_higher_valid*1.0/sum(validsgrna)
    else:
      n_higher_p=0.01
    logging.debug('higher test FDR cutoff: '+str(n_higher_p))
    # calculate gene lfc
    higher_gene_lfc=calculate_gene_lfc(args,sgrnalfc,sort_id,n_higher,sgrna2genelist,destkeys,validsgrna=validsgrna,ispos=True)
    #
    #Tracer()()
    return (n_lower_p,n_higher_p,lower_gene_lfc,higher_gene_lfc)
  else:
    return (None,None,None,None)

def rank_association_test(file,outfile,cutoff,args,adjustcutoff=True):
  if adjustcutoff: # adjust the alpha threshold to 0.05-0.5
    if cutoff<0.05:
      cutoff=0.05
    if cutoff>0.5:
      cutoff=0.5
  #rrapath='/'.join(sys.argv[0].split('/')[:-1]+["../bin/RRA"])
  rrapath='RRA'
  command=rrapath+" -i "+file+" -o "+outfile+" -p "+str(cutoff)
  if hasattr(args,'control_sgrna') and args.control_sgrna != None :
    command+=" --control "+args.control_sgrna
  if hasattr(args,'skip_gene'): 
    if args.skip_gene != None :
      for g in args.skip_gene:
        command+=" --skip-gene "+g
    else:
      command+=" --skip-gene NA --skip-gene na "
  else:
    command+=" --skip-gene NA "
  # command+=" --min-number-goodsgrna 2 "
  if hasattr(args,"additional_rra_parameters") and args.additional_rra_parameters != None:
    command+=" "+args.additional_rra_parameters+" "
  systemcall(command)


def magecktest_removetmp(prefix):
  tmpfile=[prefix+'.plow.txt',prefix+'.phigh.txt',prefix+'.gene.low.txt',prefix+'.gene.high.txt']
  for f in tmpfile:
    systemcall('rm '+f,cmsg=False)


def magecktest_parsetreatmentsfromday0(args,samplelabelindex):
  """
  Reconstruct the groups of treatment and control from --day0-label
  """
  samples=[s for s in samplelabelindex.keys()]
  day0labelstr=args.day0_label
  args.day0_label=args.day0_label.split(',')
  for dl in args.day0_label:
    if dl not in samples:
      logging.error('Label '+dl+' specified in the --day0-label option does not match count table. Please double check.')
      sys.exit(-1)
  nonday0sample=[x for x in samples if x not in args.day0_label]
  if len(nonday0sample)==0:
    logging.error('At least 1 non day0-label sample should be specified.')
    sys.exit(-1)
  args.treatment_id=nonday0sample
  args.control_id=[day0labelstr]*len(nonday0sample)


def magecktest_main(args):
  """
  Main entry for MAGeCK test function
  """
  
  warning_in_norm=False
  # stat test
  if args.subcmd == 'test':
    mapres=getcounttablefromfile(args.count_table)
  else:
    mapres=getcounttablefromfile(args.output_prefix+'.count.txt')
  cttab=mapres[0]
  sgrna2genelist=mapres[1]
  samplelabelindex=mapres[2]

  # check control sgrna

  mageckcount_checkcontrolsgrna(args,sgrna2genelist)

  if len(cttab)==0:
    sys.exit(-1)
  #nsample=len(cttab[cttab.keys()[0]])
  nsample=len(list(cttab.values())[0])

  # process day0-label
  if args.day0_label != None:
    magecktest_parsetreatmentsfromday0(args,samplelabelindex)

  # iterate control group and treatment group
  supergroup_control=args.control_id
  supergroup_treat=args.treatment_id
  # control group and treatment group labels
  labellist_control=[]
  labellist_treat=[]
  # R visualization init
  vrv=VisualRValue()
  #vrv.outprefix=args.output_prefix
  vrv.setPrefix(args.output_prefix,replace_dot=False)
  vrv.genesummaryfile=args.output_prefix+'.gene_summary.txt'
  vrv.setRmdFile()
  vrv.startRTemplate()
  vrvrnwcplabel=[]; # labels to write in rnw

  # if users specify variance estimation samples
  if args.variance_estimation_samples is not None: 
    (vargroup,vargrouplabellist)=parse_sampleids(args.variance_estimation_samples,samplelabelindex)
    logging.info('Use the following samples for variance estimation:'+','.join(vargrouplabellist))
    if len(vargroup)<2:
      logging.error('Must specify at least 2 samples for variance estimation.')
      sys.exit(-1)
  else:
    vargroup=[]
    vargrouplabellist=[]

  # loop by comparisons
  for cpindex in range(len(supergroup_treat)):
    # labels
    (treatgroup,treatgrouplabellist)=parse_sampleids(supergroup_treat[cpindex],samplelabelindex)
    treatgroup_label=str(supergroup_treat[cpindex])
    logging.info('Treatment samples:'+treatgroup_label)
    logging.info('Treatment sample index:'+','.join([str(x) for x in treatgroup]))
    labellist_treat+=[treatgroup_label]
    if supergroup_control != None:
      (controlgroup,controlgrouplabellist)=parse_sampleids(supergroup_control[cpindex],samplelabelindex)
      controlgroup_label=str(supergroup_control[cpindex]); # only for display
      logging.info('Control samples:'+controlgroup_label)
    else:
      #controlgroup=[x for x in range(nsample) if x not in treatgroup]
      #controlgrouplabellist=[samplelabelindex[x] for x in range(nsample) if x not in treatgroup]
      xls=[x for x in range(nsample) if x not in treatgroup]
      (controlgroup,controlgrouplabellist)=parse_sampleids(','.join([str(t) for t in xls]),samplelabelindex)
      controlgroup_label='rest'
      logging.info('Control samples: the rest of the samples')
    logging.info('Control sample index:'+','.join([str(x) for x in controlgroup]))
    labellist_control+=[controlgroup_label]
    #
    # convert the sample label to sample index
    if len(supergroup_treat)==1 and args.day0_label == None: 
      # directly use the prefix of the output if there's only one comparison (-t -c group), and no day0 label assigned
      cp_prefix=args.output_prefix
    else:
      # for other cases (>1 -t -c pairs, or day0label option is specified), mageck will output all individual samples vs control group results
      cp_prefix=args.output_prefix+'.'+treatgroup_label+'_vs_'+controlgroup_label # old: str(cpindex)
    #
    #
    # get the samples used
    selected_samples_id=controlgroup + treatgroup
    selected_samples_id +=[i for i in vargroup if i not in selected_samples_id]
    #
    # normalization
    cttab_sel={k:([v[i] for i in selected_samples_id]) for (k,v) in cttab.items()}; # controlgroup do not overlap with treatgroup
    norm_msg=NormMsgs()
    # perform normalization
    nttab=normalizeCounts(cttab_sel,method=args.norm_method,controlsgfile=args.control_sgrna,norm_msgs=norm_msg)
    if norm_msg.haswarning:
      warning_in_norm=True
    # write normalized counts to file
    if hasattr(args,'normcounts_to_file') and args.normcounts_to_file:
      # counts
      mageck_printdict(nttab,args,sgrna2genelist,samplelabelindex,controlgroup+treatgroup)
    # recalculate group index 
    #controlgroup_ids=list(range(len(controlgroup)))
    #treatgroup_ids=list(range(len(controlgroup),len(controlgroup+treatgroup)))
    controlgroup_ids=[i for i in range(len(selected_samples_id)) if selected_samples_id[i] in controlgroup]
    treatgroup_ids=[i for i in range(len(selected_samples_id)) if selected_samples_id[i] in treatgroup]
    vargroup_ids=[i for i in range(len(selected_samples_id)) if selected_samples_id[i] in vargroup]
    #
    #
    # for paired samples
    if hasattr(args,'paired') and args.paired:
      if len(controlgroup_ids)!=len(treatgroup_ids):
        logging.error('In paired mode, the number of control samples must be the same as the number of treatment samples.')
        sys.exit(-1)
      # treat sgRNAs separately
      new_dict={}
      new_sg_2_gene_list={}
      for (k,vi) in nttab.items():
        for ci in range(len(controlgroup_ids)):
          sgid_new=k+'_r'+str(ci)
          vi_count=[vi[controlgroup_ids[ci]],vi[treatgroup_ids[ci]]]
          vi_count+=[vi[s] for s in vargroup_ids] # if var group is specified
          new_dict[sgid_new]=vi_count
          new_sg_2_gene_list[sgid_new]=sgrna2genelist[k]
      # additional steps to replace control sgRNA ids
      if args.control_sgrna!=None:
        control_sgrna_list=[x.strip() for x in open(args.control_sgrna)]
        new_control_sgrna_list=[]
        for ci in range(len(controlgroup_ids)):
           new_control_sgrna_list+=[t+'_r'+str(ci) for t in control_sgrna_list]
        new_ctrl_file=cp_prefix+'.negctrl_ids.txt'
        nf_hd=open(new_ctrl_file,'w')
        for nli in new_control_sgrna_list:
          print(nli,file=nf_hd)
        nf_hd.close()
        args.control_sgrna=new_ctrl_file
        logging.info('The updated negative control ID file is:'+new_ctrl_file)
      # replace the original variables
      nttab=new_dict
      sgrna2genelist=new_sg_2_gene_list
      controlgroup_ids=[0]
      treatgroup_ids=[1]
      if len(vargroup_ids)>0:
        vargroup_ids=list(range(2,2+len(vargroup_ids)))
    #
    #
    # perform sgRNA test, and prepare files for gene test
    gene_as_cutoff=crispr_test(nttab, controlgroup_ids, treatgroup_ids, vargroup_ids,cp_prefix,sgrna2genelist,args)
    #
    if gene_as_cutoff[0] is not None:
      rank_association_test(cp_prefix+'.plow.txt',cp_prefix+'.gene.low.txt',gene_as_cutoff[0],args)
    if gene_as_cutoff[1] is not None:
      rank_association_test(cp_prefix+'.phigh.txt',cp_prefix+'.gene.high.txt',gene_as_cutoff[1],args,adjustcutoff=False) # update: fpr positive selection, do not adjust alpha cutoff
    # merge different files
    merge_rank_files(cp_prefix+'.gene.low.txt',cp_prefix+'.gene.high.txt',cp_prefix+'.gene_summary.txt',args,gene_as_cutoff)
    if cpindex==0: 
      if len(supergroup_treat)>1 or args.day0_label != None: # for more than 1 comparisons, set up the basic gene_summary.txt
        # for more than 1 -t -c pairs, use the results of the first pair as the default sgrna_summary.txt and gene_summary.txt 
        # if day0_label is specified, copy the files to output_prefix. This will generate files consistent with mageck-vispr
        systemcall('cp '+cp_prefix+'.gene_summary.txt '+args.output_prefix+'.gene_summary.txt',cmsg=False)
        systemcall('cp '+cp_prefix+'.sgrna_summary.txt '+args.output_prefix+'.sgrna_summary.txt',cmsg=False)
    if cpindex>0:
      if cpindex>1:
        label1=''
      else:
        if len(labellist_treat)>0:
          label1=labellist_treat[0]+'_vs_'+labellist_control[0]+'|'
        else:
          label1=''
      label2=treatgroup_label+'_vs_'+controlgroup_label+'|'
      merge_rank_summary_files(args.output_prefix+'.gene_summary.txt',cp_prefix+'.gene_summary.txt',args.output_prefix+'.gene_summary.txt',args,lowfile_prefix=label1,highfile_prefix=label2)
    # visualization: load top k genes
    # print(str(samplelabelindex))
    vrv.cplabel=treatgroup_label+'_vs_'+controlgroup_label+' neg.'
    vrvrnwcplabel+=[vrv.cplabel]
    vrv.cpindex=[2+12*cpindex+1]
    vrv.loadTopKWithExp(cp_prefix+'.gene.low.txt',nttab,sgrna2genelist,controlgrouplabellist+treatgrouplabellist)
    vrv.cplabel=treatgroup_label+'_vs_'+controlgroup_label+' pos.'
    vrvrnwcplabel+=[vrv.cplabel]
    vrv.cpindex=[2+12*cpindex+6+1]
    vrv.loadTopKWithExp(cp_prefix+'.gene.high.txt',nttab,sgrna2genelist,controlgrouplabellist+treatgrouplabellist)

    # clean the file
    if args.keep_tmp==False:
      magecktest_removetmp(cp_prefix)
      if cpindex>0:
        # systemcall('rm '+cp_prefix+'.gene_summary.txt',cmsg=False)
        # systemcall('rm '+cp_prefix+'.sgrna_summary.txt',cmsg=False) # do not remove sgrna_summary
        pass
    # end cleaning
  # end cpindex loop

  # generate pdf file
  # write to rnw file buffer
  vrv.genesummaryfile=args.output_prefix+'.gene_summary.txt'
  vrv.getGeneSummaryStat(args,isplot=False)
  vrv.comparisonlabel=vrvrnwcplabel; # replace the label field
  vrv.writeGeneSummaryStatToBuffer()
  # write to rnw and R file
  vrv.closeRTemplate()
  if hasattr(args, "pdf_report") and args.pdf_report:
    vrv.generatePDF(args.keep_tmp)
  # generate additional warning messages in terms of normalization
  if warning_in_norm:
    logging.warning('Warning messages generated during normalization. Please double check with the log file to make sure the normalization process is proper.')



