"""
MAGeCK related mathematical functions
@author:  Wei Li
@contact: li.david.wei AT gmail.com
"""
from __future__ import print_function
import sys
import math
import os
import logging


def mmedian(lst):
  """
  get the median value
  """
  sortedLst = sorted(lst)
  lstLen = len(lst)
  if lstLen==0:
    return 0.0
  index = (lstLen - 1) // 2
  
  if (lstLen % 2):
    return sortedLst[index]
  else:
    return (sortedLst[index] + sortedLst[index + 1])/2.0

def getgeomean(v):
  meanval=sum([math.log(vx+0.1,2) for vx in v])/float(len(v))
  return 2**meanval-0.1

def getMeans(matt):
  # arithmatic mean
  #meanvalue=[sum(v)/float(len(v)) for v in matt]
  # geometric mean
  # meanvalue=[getgeomean(v) for v in matt]
  # median
  meanvalue=[mmedian(v) for v in matt]
  return meanvalue

def getVars(matt):
  matt=list(matt)
  meanvalue=getMeans(matt)
  varvalue=[ sum([ (kj-meanvalue[i])*(kj-meanvalue[i])  for kj in matt[i]  ]  )/max(float(len(matt[i]))-1,1.0)    for i in range(len(meanvalue))]
  #varvalue={k:sum([ (x-meanvalue[k])*(x-meanvalue[k]) for x in v])/(float(len(v))-1)  for (k,v) in ctable.iteritems()}
  return varvalue

def leastsquare(x,y,weight=None):
  """
  least squares fitting
  coefficients from y= a+bx
  return (b,a)
  reference: http://mathworld.wolfram.com/LeastSquaresFitting.html
  For weighted least square: http://goo.gl/pGpTZ6
  """
  n=len(x)
  if n != len(y):
    logging.error('Unequal length of vectors of x and y in least square')
    sys.exit(-1)
  if weight is None:
    sy=sum(y)
    sx=sum(x)
    sx2=sum([t*t for t in x])
    sxy=sum([x[i]*y[i] for i in range(n)])
    a=(sy*sx2-sx*sxy)/(n*sx2-sx*sx)
    b=(n*sxy-sx*sy)/(n*sx2-sx*sx)
    return (b,a)
  else:
    nw=sum(weight)
    sy=sum([y[i]*weight[i] for i in range(n)])
    sx=sum([x[i]*weight[i] for i in range(n)])
    sx2=sum([x[i]*x[i]*weight[i] for i in range(n)])
    sxy=sum([x[i]*y[i]*weight[i] for i in range(n)])
    a=(sy*sx2-sx*sxy)/(nw*sx2-sx*sx)
    b=(nw*sxy-sx*sy)/(nw*sx2-sx*sx)
    return (b,a)


def getnormcdf(x,lowertail=True):
  """
  Get the normal CDF function. used to calculate p-value
  """
  # ax=math.fabs(x)
  #axv=math.erfc(x/(2**0.5))/2; # higher tail
  if lowertail==False:
    #return axv
    return math.erfc(x/(2**0.5))/2
  else:
    #return 1-axv
    return math.erfc(-x/(2**0.5))/2
  #if (x>0 and lowertail==False) or (x<0 and lowertail==True):
  #  return axv
  #else:
  #  return 1-axv

def getNormalPValue(mean0,var0,mean1, lower=False):
  """
  Use truncated normal distribution to calculate the pvalue
  """
  # use ttmean to calculate the pvalue
  n=len(mean0)
  minmean1=min([x for x in mean1 if x>0])
  mean1_adj=[(lambda x: x if x >minmean1 else minmean1)(t) for t in mean1]
  # first, convert to standard normal distribution values
  t_theta=[(mean1_adj[i]-mean0[i])/math.sqrt(var0[i]) for i in range(n)]
  t_theta_0=[(0.0-mean0[i])/math.sqrt(var0[i]) for i in range(n)]
  #
  t_p=[getnormcdf(x,lowertail=lower) for x in t_theta]
  t_p_0=[getnormcdf(x,lowertail=True) for x in t_theta_0]
  if lower==True:
    return [(t_p[i]-t_p_0[i])/(1-t_p_0[i]) for i in range(n)]
  else:
    return [t_p[i]/(1-t_p_0[i]) for i in range(n)]


def getNBPValue(mean0,var0,mean1, lower=False,log=False):
  """
  Use negative binomial to calculate p-value
  Reference:
  http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html#scipy.stats.nbinom
  """
  from scipy.stats import nbinom
  n=len(mean0)
  nb_p=[mean0[i]/var0[i] for i in range(n)]; # consisitent with R
  nb_n0=[mean0[i]*mean0[i]/(var0[i]-mean0[i]) for i in range(n)]
  nb_n=[ (lambda t: t if t>=1 else 1)(x) for x in nb_n0]
  #
  if lower==True:
    if log==False:
      nb_p_low=nbinom.cdf(mean1,nb_n,nb_p)
    else:
      nb_p_low=nbinom.logcdf(mean1,nb_n,nb_p)
    return list(nb_p_low)
  else:
    if log==False:
      nb_p_low=nbinom.sf(mean1,nb_n,nb_p)
    else:
      nb_p_low=nbinom.logsf(mean1,nb_n,nb_p)
    return list(nb_p_low)


# an accurate estimation of qnorm 

class QNormConverter:
  """
  A simple converter to calculate qnorm (or scipy.stats.norm.ppf) based on pre_computed value
  """
  logp_vec=[]
  norm_qnorm_vec=[]
  qmin=None
  qmax=None
  qstep=None
  def  __init__(self):
    n=0
    qnormfile=os.path.join(os.path.dirname(__file__),'qnorm_table.txt')
    try:
      for lines in open(qnormfile,'r'):
        n+=1
        if n==1:
          continue
        field=lines.strip().split()
        self.logp_vec+=[float(field[0])]
        self.norm_qnorm_vec+=[float(field[1])]
    except IOError:
      logging.error('Cannot open qnorm_table.txt in the default mageck source folder: '+qnormfile+' for QC. Please double check your file location and format.')
      sys.exit(-1)
    self.qmin=min(self.logp_vec)
    self.qmax=max(self.logp_vec)
    self.qstep=self.logp_vec[1]-self.logp_vec[0]
  
  def get_qnorm(self,p,islog=False):
    """
    Look up convert table
    """
    if type(p) is not list:
      p_list=[p]
    else:
      p_list=p
    if islog:
      lp=p_list
    else:
      lp=[math.log(x) for x in p_list]
    # calculate index
    vec_len=len(self.logp_vec)
    ind_min=[int(math.floor((x-self.qmin)/self.qstep)) for x in lp]
    ind_min=[min(max(x,0),vec_len-2) for x in ind_min]
    frac_value=[ (lp[i]-self.logp_vec[ind_min[i]])/self.qstep for i in range(len(lp))]
    frac_value=[min(max(x,0.0),1.0) for x in frac_value]
    interp_value=[ self.norm_qnorm_vec[ind_min[i]]*(1.0-frac_value[i]) + self.norm_qnorm_vec[ind_min[i]+1]*frac_value[i] for i in range(len(lp))]
    interp_value_ret=[interp_value[i] if (lp[i] >= self.qmin and lp[i]<=self.qmax) else None for i in range(len(lp))]
    if type(p) is not list:
      return interp_value_ret[0]
    else:
      return interp_value_ret
    
    
