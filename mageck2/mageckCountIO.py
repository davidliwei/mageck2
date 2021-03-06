""" Processing various file types for MAGeCK count
"""

from __future__ import print_function

import sys
import math
import logging
import string


def mageckcount_gini(x):
  '''
  Return the Gini index of an array
  Calculation is based on http://en.wikipedia.org/wiki/Gini_coefficient
  '''
  xs=sorted(x)
  n=len(xs)
  gssum=sum([ (i+1.0)*xs[i] for i in range(n)])
  ysum=sum(xs)
  if ysum==0.0:
    ysum=1.0
  gs=1.0-2.0*(n-gssum/ysum)/(n-1)
  return gs

def mageckcount_search_variable_region(seqlist):
  '''
  Search variable region of the sequence, possibly containing UMIs or the second pair
  Parameter:
  -----------
  seqlist
    A list of sequence
  Return value
  ----------
  (start,end)
    A variable region containing guides/UMIs in the seqlist
  '''
  count_freq={}
  for seq in seqlist:
    sequ=seq.upper()
    for i in range(len(seq)):
      if i not in count_freq:
        count_freq[i]=[0.0,0.0,0.0,0.0]  # A T G C
      if sequ[i] == 'A':
        count_freq[i][0]+=1
      if sequ[i] == 'T':
        count_freq[i][1]+=1
      if sequ[i] == 'G':
        count_freq[i][2]+=1
      if sequ[i] == 'C':
        count_freq[i][3]+=1
  var_start=-1
  var_end=-1
  #print(count_freq)
  logging.info('--count-freq data for UMI search (A/T/G/C):')
  for x in range(max(count_freq.keys())):
    logging.info('  '+str(x)+'\t'+'/'.join([str(s) for s in count_freq[x]]))
  
  # search for start
  for x in range(max(count_freq.keys())):
    freq_vec=count_freq[x]
    freq_t_sum=max(sum(freq_vec),1)
    freq_f=[y/freq_t_sum for y in freq_vec]
    if freq_t_sum*1.0/len(seqlist) > 0.9 and max(freq_f)>0 and max(freq_f)<0.6:
      var_start=x
      break
  if var_start>=0:
    for x in range(var_start+1,max(count_freq.keys())):
      freq_vec=count_freq[x]
      freq_t_sum=max(sum(freq_vec),1)
      freq_f=[y/freq_t_sum for y in freq_vec]
      if freq_t_sum*1.0/len(seqlist) > 0.9 and max(freq_f)>0 and max(freq_f)<0.6:
        pass
      else:
        var_end=x
        break
  return (var_start,var_end)
  


def mageckcount_trim5_auto(filename,args,genedict, revcomp=False,is_second_pair=False,no_search=False):
  '''
  Automatically determine the trim 5 length in one fastq file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  genedict
    Library
  revcomp
    Shall the reads be reverse complemented?
  is_second_pair
    Is this search on the second read of the paired read?
  no_search
    Do not search for guide position (return the entire sequence in remainingseq_list). Note that revcomp won't happen

  Return value
  -----------
  candidatetrim5
    A list of possible trim-5 value
  remainingseq_list
    A list of remaining sequences after sgrna

  '''
  # ctab={}
  logging.info('Determining the trim-5 length of FASTQ file '+filename+'...')
  if len(genedict)==0:
    logging.error('Library file must be provided to automatically determine the trim-5 length in fastq file.')
  # checking possible sgRNA length
  lengthpool={}
  for k in genedict.keys():
    if len(k) not in lengthpool:
      lengthpool[len(k)]=0
    lengthpool[len(k)]+=1
  lengthpoolkeys=sorted(lengthpool.keys(),reverse=True)
  minlengthpool=min(lengthpool)
  trimkey={}
  readlenkey={}
  logging.info('Possible gRNA lengths:'+','.join([str(t) for t in lengthpoolkeys]))
  if filename.upper().endswith('.GZ'):
    import gzip
    openobj=gzip.open(filename,'rt')
  else:
    openobj=open(filename)
  # determine whether trim-5 value is provided in arguments, or will need to determine here
  trim_5_provided=None
  if is_second_pair:
    trim_5_provided=None # trim-5 only works for first pair
  elif args.trim_5.upper() == 'AUTO' or args.trim_5.upper() == 'AUTOTEST':
    trim_5_provided=None
  else:
    try:
      trim_5_provided=[int(x) for x in args.trim_5.split(',')]
      logging.info('Specified trimming length:'+','.join([str(x) for x in candidate_trim5]))
    except ValueError:
      logging.error('Integer values must be specified in --trim-5 option')

  
  nline=0
  maxline=100000 # maximum reads tested
  nreadcount=0
  remainingseq_list=[]
  for line in openobj:
    # line=line.encode('utf-8')
    nline=nline+1
    if nline%4 == 2:
      nreadcount+=1
      if nreadcount>maxline:
        break
      if nreadcount%100000==1:
        logging.info('Processing '+str(round(nreadcount/1000000))+ 'M reads ...')
      fseq=line.strip()
      if no_search:
        remainingseq_list+=[fseq]
        continue
      if len(fseq) not in readlenkey:
        readlenkey[len(fseq)]=0
      readlenkey[len(fseq)]+=1
      # check length
      # for l in lengthpool.keys():
      trim_test_range=[tl_r for tl_r in range(len(fseq)-minlengthpool+1)]
      if trim_5_provided is not None:
        trim_test_range=trim_5_provided
      for triml in trim_test_range:
        fseqtrim=fseq[triml:]
        findrecord=False
        for l in lengthpoolkeys: # iterate all possible lengths
          testl=l
          if len(fseqtrim)<testl:
            continue
          fseqc=fseqtrim[:testl]
          if fseqc.count('N')>0 and args.count_n==False:
            continue
          if revcomp:
            fseqc=mageckcount_revcomp(fseqc)
          if fseqc not in genedict:
            continue
          else:
            # find the record
            if triml not in trimkey:
              trimkey[triml]=0
            trimkey[triml]+=1
            findrecord=True
            remainingseq=fseqtrim[testl:] # add the remaining sequences
            remainingseq_list+=[remainingseq]
            break
        # end for l
        if findrecord:
          break
      # end for triml 
  # end for line 
  openobj.close()

  if no_search:
    return ([],remainingseq_list)
  
  keysorted=sorted(trimkey.items(),key=lambda z: z[1], reverse=True)
  totalmappedreads=sum([x[1] for x in trimkey.items()])
  totalfrac=totalmappedreads*1.0/nreadcount
  logging.info('Read length:'+','.join([str(x) for x in readlenkey.keys()]))
  logging.info('Total tested reads: '+str(nreadcount)+', mapped: '+str(totalmappedreads)+ '('+str(totalfrac)+')')
  if totalfrac < 0.001:
    logging.error('Cannot automatically determine the --trim-5 length. Only ' + str(totalfrac*100 )+ ' % of the reads can be identified. ')
    sys.exit(-1)
  elif totalfrac < 0.5:
    logging.warning('Only ' + str(totalfrac*100 )+ ' % of the reads can be identified. ')
  candidatetrim5=[]
  candidatetrim5frac=[]
  lastfrac=0.01
  logging.info('--trim-5 test data: (trim_length reads fraction)')
  for ks in keysorted:
    ksfrac=ks[1]*1.0/totalmappedreads
    logging.info('\t'.join([str(x) for x in [ks[0],ks[1],ksfrac]]))
  for ks in keysorted:
    ksfrac=ks[1]*1.0/totalmappedreads
    # logging.info('\t'.join([str(x) for x in [ks[0],ks[1],ksfrac]]))
    if lastfrac>ksfrac*3.0:
      break
    candidatetrim5+=[ks[0]]
    candidatetrim5frac+=[ksfrac]
    if sum(candidatetrim5frac)>0.75:
      break
    if len(candidatetrim5)>0 and ksfrac<0.05:
      break
    lastfrac=ksfrac
  logging.info('Auto determination of trim5 results: '+','.join([str(x) for x in candidatetrim5]))
  return (candidatetrim5,remainingseq_list)
  # return 0



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

def mageckcount_search_trim_and_sglen(args,fseq0,genedict,ctab,ctab_umi,candidate_trim5,lengthpoolkeys,revcomp=False):
  '''
  within one line, search for best matches of 5' trimming length
  Parameters
  ------------
  args
    arguments
  fseq0
    sequence
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  ctab
    A dictionary of sgRNA sequence and count
  ctab_umi
    a {guide:{UMI:count}} dictionary structure
  candidate_trim5
    A list of possible trim5 values
  lengthpoolkeys
    possible sgrna lengths
  revcomp
    Whether to reverse complement the sequence?
    
  Return value:
  ------------
    findrecord
      True if a record is found, false otherwise
    guideseq
      The guide sequence found; none if not found
  '''
  findrecord=False
  fseq=''
  matched_seq=None
  for triml in candidate_trim5:
    # check length
    fseq=fseq0[triml:]
    # for l in lengthpool.keys():
    if len(genedict)==0:
      # WL: looks like this branch will consider the rest of the reads as guides 
      if len(fseq)<args.sgrna_len:
        continue
        #fseq=fseq[:args.sgrna_len]
      if revcomp:
         fseq = mageckcount_revcomp(fseq)
      if (fseq.count('N')>0 ) and args.count_n==False:
        continue
      if fseq not in ctab:
        ctab[fseq]=0
      ctab[fseq]=ctab[fseq]+1
    else:
      for l in lengthpoolkeys: # iterate all possible lengths
        testl=l
        if len(fseq) < testl:
          continue
        fseqc = fseq[:testl]
        if revcomp:
          fseqc = mageckcount_revcomp(fseqc)
        # Only use single end fastq file
        if (fseqc.count('N')>0 ) and args.count_n == False :
          continue
        if fseqc not in genedict:
          continue
        # Count the number of sgRNA
        if fseqc not in ctab:
          ctab[fseqc]=0
        ctab[fseqc]=ctab[fseqc]+1
        findrecord=True
        matched_seq=fseqc
        if args.umi=='firstpair':
          # extract UMIs from the first pair
          if args.umi_start == -1 or args.umi_end == -1:
            logging.error('Error: umi-start and umi-end has to be greater than zero')
          umiseq = fseq[(testl+args.umi_start):(testl+args.umi_end)]
          if fseqc not in ctab_umi:
            ctab_umi[fseqc]={}
          if umiseq not in ctab_umi[fseqc]:
            ctab_umi[fseqc][umiseq]=0
          ctab_umi[fseqc][umiseq]=ctab_umi[fseqc][umiseq]+1
        # nmappedcount+=1
        break
      # end for l
      if findrecord:
        break
    if findrecord:
      break
  # end for loop of trim-5
  # save unmapped file
  if args.unmapped_to_file and findrecord==False:
    if len(fseq)>=args.sgrna_len:
      fseqc=fseq[:args.sgrna_len]
      # fseqd=mageckcount_revcomp(fseqr[:args.sgrna_len])
      if fseqc.count('N') > 0 and args.count_n == False :
        pass
      else:
        if fseqc not in ctab:
          ctab[fseqc]=0
        ctab[fseqc]=ctab[fseqc]+1
    
  return (findrecord, matched_seq)


def mageckcount_determine_umi_location(args,pairedfile,genedict,remainingseq_list):
  '''
  Determine the parameters for UMI and paired guide
  Parameters
  -----------
  Return value
  -----------
  '''
  (umi_start,umi_end)=mageckcount_search_variable_region(remainingseq_list)
  if umi_start >= 0 and umi_end >= 0:
    logging.info('UMI found in the first read. Position (after guide): '+str(umi_start)+'-'+str(umi_end))
    args.umi_start=umi_start
    args.umi_end=umi_end
    args.umi='firstpair'
  else:
    if pairedfile != None:
      logging.info('Search for UMI in the second read...')
      (candidate_trim5_paired, remainingseq_list_pair) = mageckcount_trim5_auto(pairedfile, args, genedict, revcomp=True,is_second_pair=True, no_search=True)
      (umi_start_2,umi_end_2)=mageckcount_search_variable_region(remainingseq_list_pair)
      if umi_start_2 >= 0 or umi_end_2 >= 0:
        logging.info('UMI found in the second read. Position: '+str(umi_start_2)+'-'+str(umi_end_2))
        args.umi_start_2=umi_start_2
        args.umi_end_2=umi_end_2
        args.umi='secondpair'
      else:
        logging.error('Search failed for UMI in the second read.')
        sys.exit(-1)
    else:
     logging.error('Search failed for UMI in the first read.')
     sys.exit(-1)
  return 0


def mageckcount_processonefile(filename,args,ctab,ctab_umi,genedict,datastat,pairedfile, adjust):
  '''
  Go through one fastq file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A {sequence:count} dictionary
  ctab_umi
    A {sequence:{umi:count}} dictionary
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})
  pairedfile
    Paired fastq filename to be sequence
  adjust
    Adjust for paired counts

  Return value
  -----------
  datastat
    a dictionary structure of statistics
  '''
  # ctab={}
  nline=0
  logging.info('Parsing FASTQ file '+filename+'...')
  nreadcount=0
  nmappedcount=0
  # checking possible trimming lengths
  candidate_trim5=[0]
  remainingseq_list=[]
  candidate_trim5_paired=[0]
  remainingseq_list_pair=[]
  (candidate_trim5,remainingseq_list)=mageckcount_trim5_auto(filename,args,genedict)
  # parameters for UMI search

  revised_umi_auto=False # a marker to retore args.umi later
  #ctab_umi={} # a {guide:{UMI:count}} dictionary structure
  if args.umi != 'none':
    if args.umi=='auto':
      # automatically search for UMIs
      revised_umi_auto=True# a marker to retore args.umi later
      mageckcount_determine_umi_location(args,pairedfile,genedict,remainingseq_list)
  elif args.pairguide != 'none':
    if args.pairguide=='auto':
      mageckcount_determine_umi_location(args,pairedfile,genedict,remainingseq_list)
      # copy the configurations from umi to pg, and retore the configurations of umi
      if args.umi=='firstpair':
        args.pg_start=args.umi_start
        args.pg_end=args.umi_end
        args.pairguide='firstpair'
        args.umi='none'
        args.umi_start=-1
        args.umi_end=-1
      elif args.umi=='secondpair':
        args.pg_start_2=args.umi_start_2
        args.pg_end_2=args.umi_end_2
        args.pairguide='secondpair'
        args.umi='none'
        args.umi_start_2=-1
        args.umi_end_2=-1
  else:
    if pairedfile != None:
      (candidate_trim5_paired, remainingseq_list_pair) = mageckcount_trim5_auto(pairedfile, args, genedict, revcomp=True,is_second_pair=True)

  # do not allow multiple trim_5 options without library file
  if len(candidate_trim5)>1 and len(genedict)==0:
    logging.error('A library file has to be provided if multiple trimming lengths are specified in --trim-5 option.')
    sys.exit(-1)
  if len(candidate_trim5)>1 and args.unmapped_to_file:
    logging.error('The --unmapped-to-file option is not allowed if multiple trimming lengths are specified in --trim-5 option.')
    sys.exit(-1)
  # checking possible sgRNA length
  lengthpool={}
  for k in genedict.keys():
    if len(k) not in lengthpool:
      lengthpool[len(k)]=0
    lengthpool[len(k)]+=1
  lengthpoolkeys=sorted(lengthpool.keys(),reverse=True)
  logging.info('Possible gRNA lengths:'+','.join([str(t) for t in lengthpoolkeys]))
  # Open fastq file
  if filename.upper().endswith('.GZ'):
    import gzip
    openobj=gzip.open(filename,'rt')
    if pairedfile != None:
      openpaired = gzip.open(pairedfile, 'rt')
  else:
    openobj=open(filename)
    if pairedfile != None:
      openpaired = open(pairedfile)

  # Begin to read fastq file
  for line in openobj:
    # line=line.encode('utf-8')
    pairline=line
    if pairedfile != None:
      pairline = openpaired.readline()
    nline=nline+1
    if nline%4 != 2:
      continue
    nreadcount+=1
    if nreadcount%1000000==1:
      logging.info('Processing '+str(round(nreadcount/1000000))+ 'M reads ..')
      if nreadcount>1000000 and hasattr(args, 'test_run') and args.test_run:
        break
    fseq0=line.strip()
    fseq1=pairline.strip()
    #if args.trim_5 >0:
    #  fseq=fseq0[args.trim_5:]
    # search the firstpair
    (firstpair_found, matched_seq)=mageckcount_search_trim_and_sglen(args,fseq0,genedict,ctab,ctab_umi,candidate_trim5,lengthpoolkeys)
    if firstpair_found == True:
      nmappedcount+=1
    if args.umi!='none':
      # now search for UMIs on the second pair
      if args.umi=='secondpair':
        # extract UMIs from the second pair
        if args.umi_start_2 == -1 or args.umi_end_2 == -1:
          logging.error('Error: umi-start and umi-end has to be greater than zero')
        umiseq = fseq1[args.umi_start_2:args.umi_end_2]
        if matched_seq not in ctab_umi:
          ctab_umi[matched_seq]={}
        if umiseq not in ctab_umi[matched_seq]:
          ctab_umi[matched_seq][umiseq]=0
        ctab_umi[matched_seq][umiseq]=ctab_umi[matched_seq][umiseq]+1
    else:
      if firstpair_found == False  and adjust==True and pairedfile != None:
        # if not found in the first pair, search the second pair
        (secondpair_found,matched_seq)=mageckcount_search_trim_and_sglen(args,fseq1,genedict,ctab,ctab_umi,candidate_trim5_paired,lengthpoolkeys,revcomp=True)
        if secondpair_found == True:
          nmappedcount+=1

  # 
  logging.info('Total: '+str((nreadcount))+ '.')
  logging.info('Mapped: '+str((nmappedcount))+ '.')
  openobj.close()
  if pairedfile != None:
    openpaired.close()
  # calculate statistics
  datastat['reads']=nreadcount
  # check if a library is provided
  if len(genedict)==0:
    datastat['mappedreads']=0
    datastat['totalsgrnas']=0
    datastat['zerosgrnas']=0
    datastat['giniindex']=1
  else:
    nmapped=0
    nrdcnt=[]
    for (k,v) in ctab.items():
      if k in genedict:
        nmapped+=v
        nrdcnt+=[math.log(v+1.0)]
    nzerosg=0
    for (k,v) in genedict.items():
      if k not in ctab:
        nzerosg+=1
        nrdcnt+=[math.log(0.0+1.0)]
    # logging.info('mapped:'+str(nmapped))
    datastat['mappedreads']=nmapped
    datastat['totalsgrnas']=len(genedict);
    datastat['zerosgrnas']=nzerosg
    datastat['giniindex']=mageckcount_gini(nrdcnt)
    #
    #print('Mapped: '+str((nmappedcount))+ '...')
    #print('Mapped 2: '+str((nmapped))+ '...')
  #return ctab
  if False and revised_umi_auto:# a marker to retore args.umi the default values later; currently off
    args.umi='auto'
    args.umi_start=-1
    args.umi_end=-1
    args.umi_start_2=-1
    args.umi_end_2=-1
  return 0




def mageckcount_processonefile_bam(filename,args,ctab,genedict,datastat):
  '''
  Go through bam file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A dictionary of sgRNA sequence and count
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})

  Return value
  -----------
  datastat
    a dictionary structure of statistics
  
  # important note for the alignment
  1. Make sure 5' and 3' adapters are properly removed before mapping. Either use cutadapt or --trim5/--trim3 option in bowtie2.
  2. Make sure no reverse-complement mapping is allowed; otherwise, there will be multiple mappings for sgRNAs whose sequnces are reverse complemented. In bowtie2, --no-rc parameter should be specified.
  3. Carefully check the alignment strategy used in the aligner. Some sequences may map to multiple sgRNAs and are counted multiple times.
  4. When building index, some aligners (like bowtie2) will remove sgRNAs with identical sequence. This will create some warning messages "sgRNA in the BAM file does not match th provided library file".

  Reference: 
  BAM specification 
  https://samtools.github.io/hts-specs/SAMv1.pdf 
  Reference: Tao Liu's MACS2 
  https://github.com/taoliu/MACS/blob/master/MACS2/IO/Parser.pyx
  '''
  import sys 
  import gzip 
  import io 
  import math 
  import struct 
  from struct import unpack 
  """ 
  Encode table  
  ACMGRSVTWYHKDBNN -> [0,15] 
  Code int bit bit(reverse) int(reverse) 
    0  0000    
  A 1  0001    
  C 2  0010    
  M 3  0011    
  G 4  0100    
  R 5  0101    
  S 6  0110    
  V 7  0111    
  T 8  1000    
  W 9  1001    
  Y 10 1010    
  H 11 1011    
  K 12 1100    
  D 13 1101    
  B 14 1110    
  N 15 1111    
  """ 
  encodetable='NACMGRSVTWYHKDBN' 
  nline=0
  logging.info('Parsing BAM file '+filename+'...')
  nreadcount=0
  
  # start processing the bam file
  # open bam file 
  fhd=io.BufferedReader( gzip.open( filename, mode='rb' ) ) 
  # check the first 3 character must be BAM 
  fhd.seek(0) 
  magic_header = fhd.read( 3 ) 
  if magic_header.decode("utf-8")!= "BAM": 
    logging.error('Error: not recognized BAM file: '+filename+', header:'+magic_header.decode("utf-8")) 
    sys.exit(-1) 
  # check header 
  fhd.seek( 4 ) 
  header_len =  unpack( '<i', fhd.read( 4 ) )[ 0 ] 
  fhd.seek( header_len + fhd.tell() ) 
  # next, get chromosome, and check whether it matches the given genedict
  genedict_sgid={v[0]:k for (k,v) in genedict.items()}
  nc = unpack( '<i', fhd.read( 4 ) )[ 0 ] 
  refnames=['']*nc 
  refnameslen=[0]*nc 
  nwarningsg=0
  for x in range( nc ): 
    # read each chromosome name 
    nlength = unpack( '<i' , fhd.read( 4 ) )[ 0 ] 
    refstr=fhd.read( nlength ).decode("utf-8") 
    # jump over chromosome size, we don't need it 
    refstrlen = unpack( '<i', fhd.read( 4 ) )[ 0 ] 
    #fhd.seek( fhd.tell() + 4 ) 
    #print(refstr+':'+str(refstrlen)) 
    refstr=refstr[:-1]
    refnames[x]=refstr 
    refnameslen[x]=refstrlen 
    if refstr not in genedict_sgid:
      # logging.warning('sgRNA ID '+ refstr+' in the BAM file does not not match the provided library file. Please double check.')
      if nwarningsg<3:
        logging.warning('sgRNA ID '+ refstr+' in the BAM file is missing in the provided library file (it may have duplicated sequences with other sgRNAs). Please double check.')
      nwarningsg+=1
  logging.info(str(nc)+' references detected in the BAM file.') 
  if nwarningsg>3:
    logging.warning('Total sgRNAs that are missing in the library (may be due to duplicated sequences):'+str(nwarningsg))
  nwarningsg=0
  #
  # next, iterate the bam file
  while True: 
    nline=nline+1
    if nline%1000000==1:
      logging.info('Processing '+str(round(nline/1000000))+ 'M records ..')
      if nline>1000000 and hasattr(args, 'test_run') and args.test_run:
        break
    # 

    tmpdata=fhd.read(4) 
    if len(tmpdata)==0: 
      break 
    entrylength = unpack( '<i', tmpdata )[ 0 ] 
    data = fhd.read( entrylength )
    # Skip paired reads when the first read is mapped.
    next_pos = unpack('<i', data[24:28])[0]
    # print(str(nreadcount)+":"+str(next_pos))
    # if next_pos<(nreadcount+1):
    if next_pos >0 and next_pos<(nreadcount+1):
      continue
    nreadcount += 1
    # refid, position 
    refid=unpack( '<i', data[:4] )[ 0 ] 
    if refid == -1: 
      # didn't find any match
      refstr='*' 
    else: 
      # find matches
      refstr=refnames[refid] 
      if refstr not in genedict_sgid:
        # logging.warning('sgRNA ID: '+refstr+' is not present in the library file. Please double-check the consistency between library file and SAM/BAM file.')
        if nwarningsg<3:
          logging.warning('sgRNA ID: '+refstr+' is not present in the library file. Please double-check the consistency between library file and SAM/BAM file.')
        nwarningsg+=1
        continue
      fseqc=genedict_sgid[refstr]
      if fseqc not in ctab:
        ctab[fseqc]=0
      ctab[fseqc]=ctab[fseqc]+1

    # other fields in BAM file; not used
    if False:
      position=unpack( '<i', data[4:8] )[ 0 ] 
      bin_mq_nl=unpack( '<i', data[8:12] )[ 0 ] 
      read_name_len=bin_mq_nl&0x000000ff  
      # print('name length:'+str(read_name_len)) 
      flag_nc=unpack( '<i', data[12:16] )[ 0 ] 
      n_cigar_op=flag_nc&0x0000ffff 
      flag=flag_nc>>16 
      # length 
      readlen = unpack( '<i', data[16:20] )[ 0 ] 
      next_refID=unpack( '<i', data[20:24] )[ 0 ] 
      next_pos=unpack( '<i', data[24:28] )[ 0 ] 
      tlen=unpack( '<i', data[28:32] )[ 0 ] 
      read_name=data[32:32+read_name_len] 
      # cigar 
      tstart=32+read_name_len 
      tend=tstart+n_cigar_op*4 
      # the following is the only 1st int of cigar string 
      cigarstr=unpack('<i',data[tstart:tstart+4]) 
      # seq, encoded 
      tstart=tend 
      tend=tstart+int(math.ceil((readlen+1)/2)) 
      seq=data[tstart:tend] 
      # quality 
      tstart=tend 
      tend=tstart+readlen 
      phredqual=data[tstart:tend] 
      if nline<10 and False: # only for debug purposes 
        logging.info('refid: '+str(refid)+' '+refstr+', position:'+str(position)+', name:'+read_name) 
        # decode the sequence 
        nseqcode=int(math.floor((readlen+1)/2)) 
        seqstring='' 
        for i in range(nseqcode): 
          iit=seq[i] 
          iit_int=unpack('<B',iit)[0] 
          k1=iit_int>>4 
          k2=iit_int&0xf 
          seqstring+=encodetable[k1] 
          if readlen %2==1 and i==nseqcode-1: 
            pass 
          else: 
            seqstring+=encodetable[k2] 
        # end for
        logging.info(seqstring) 
      # end if
    # end if FALSE
  # end while
  if nwarningsg>3:
    logging.warning('Total records that are not in the library:'+str(nwarningsg))
  fhd.close() 
  
  # calculate statistics
  datastat['reads']=nreadcount
  nmapped=0
  nrdcnt=[]
  for (k,v) in ctab.items():
    if k in genedict:
      nmapped+=v
      nrdcnt+=[math.log(v+1.0)]
  nzerosg=0
  for (k,v) in genedict.items():
    if k not in ctab:
      nzerosg+=1
      nrdcnt+=[math.log(0.0+1.0)]
  logging.info('mapped:'+str(nmapped))
  datastat['mappedreads']=nmapped
  datastat['totalsgrnas']=len(genedict);
  datastat['zerosgrnas']=nzerosg
  datastat['giniindex']=mageckcount_gini(nrdcnt)
  #return ctab
  return 0


def mageckcount_processonefile_sam(filename,args,ctab,genedict,datastat):
  '''
  Go through sam file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A dictionary of sgRNA sequence and count
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})

  Return value
  -----------
  datastat
    a dictionary structure of statistics
  
  # important note for the alignment
  # Please see the notes for bam files

  Reference: 
  BAM specification 
  https://samtools.github.io/hts-specs/SAMv1.pdf 
  Reference: Tao Liu's MACS2 
  https://github.com/taoliu/MACS/blob/master/MACS2/IO/Parser.pyx
  '''
  import sys 
  import gzip 
  import io 
  import math 
  import struct 
 
  nline=0
  logging.info('Parsing SAM file '+filename+'...')
  nreadcount=0
  
  # the {sgRNA_id: sequence} directory
  genedict_sgid={v[0]:k for (k,v) in genedict.items()}
  # whether the warning of particular sgRNA id is already present
  sgrnaid_haswarning={}
  # start processing the sam file
  for line in open(filename):
    if line[0]=='@':
      continue
    nline=nline+1
    if nline%1000000==1:
      logging.info('Processing '+str(round(nline/1000000))+ 'M records ..')
      if nline>1000000 and hasattr(args, 'test_run') and args.test_run:
        break
    #
    field=line.strip().split()
    # refid, position 
    refid=field[2]
    # Pnext=(field[7])
    Rnext=(field[8])
    Pnext=int(field[7])
    POSv=int(field[3])
    # Skip paired reads when the first read is mapped.
    if Pnext > 0 and Pnext<(nreadcount+1):
    ### if Rnext!="*" and Pnext != 0 and Pnext<(POSv):
      continue
    nreadcount += 1

    if refid != "*": 
      # find matches
      refstr= refid
      if refstr not in genedict_sgid:
        if refstr not in sgrnaid_haswarning:
          # logging.warning('sgRNA ID: '+refstr+' in the SAM file does not match the provided library file. Please double-check.')
          sgrnaid_haswarning[refstr]=1
        continue
      fseqc=genedict_sgid[refstr]
      if fseqc not in ctab:
        ctab[fseqc]=0
      ctab[fseqc]=ctab[fseqc]+1
  # end while
  # calculate statistics
  datastat['reads']=nreadcount
  nmapped=0
  nrdcnt=[]
  for (k,v) in ctab.items():
    if k in genedict:
      nmapped+=v
      nrdcnt+=[math.log(v+1.0)]
  nzerosg=0
  for (k,v) in genedict.items():
    if k not in ctab:
      nzerosg+=1
      nrdcnt+=[math.log(0.0+1.0)]
  logging.info('mapped:'+str(nmapped))
  logging.warning('sgRNAs not in the library (may be due to duplicated sequences):'+str(len(sgrnaid_haswarning)))
  datastat['mappedreads']=nmapped
  datastat['totalsgrnas']=len(genedict);
  datastat['zerosgrnas']=nzerosg
  datastat['giniindex']=mageckcount_gini(nrdcnt)
  #return ctab
  return 0



