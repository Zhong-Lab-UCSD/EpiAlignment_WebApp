# This version is for EpiAlignment Web APP.

import sys, argparse
from xplib.Annotation import Bed
from xplib import DBI
from subprocess import Popen, PIPE
import string
from multiprocessing import *
from time import time

rev_table=string.maketrans('ACGTacgtN', 'TGCATGCAN')

def ParseArg():
  p=argparse.ArgumentParser(description="Generate fastq file for EpiAlignment" )
  p.add_argument("q_region",type=str, nargs='+', help="coordinates of query regions. chr:start-end:strand. Note that the order should be the same as species.")
  p.add_argument("-s", "--species", nargs=2, type=str, default=["hg38","mm10"], help="A list of genome assembly. Note that the list should have the same order as the bedgraph files")
  p.add_argument("-f","--fasta_path", type=str, default="Annotation/Sequence", help="Path to the genome sequence files.")
  p.add_argument("--histone", nargs="+",type=str, default=["H3K4me3"], help="Name of histone modifications. Note that the list should have the same order as histone.")
  p.add_argument("--bg", type=str, help="file name. The file contains paths to ChIP-Seq peak calling files.")
  p.add_argument("--s_path",type=str, default="samtools", help="path of samtools")
  p.add_argument("-p","--p_num",type=int, default=5, help="Number of processes.")
  p.add_argument("-o","--output",type=str,help="output file name.")
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args() 

def overlap(bed1,bed2):
    if (bed1.stop > bed2.start) and (bed1.start < bed2.stop):
        return min(bed1.stop,bed2.stop) - max(bed1.start,bed2.start) + 1
    else:
        return False

def ReadHistones(fhist_name):
  sp1 = []
  sp2 = []
  with open(fhist_name, "r") as fhist:
    fhist.readline()
    while True:
      line = fhist.readline().strip()
      if "@" in line:
        break
      line = line.split()
      sp1.append(DBI.init(line[0], "bed"))
    while True:
      line = fhist.readline().strip()
      if line == "":
        break
      line = line.split()
      sp2.append(DBI.init(line[0], "bed"))
 
  return sp1, sp2


def revcomp(seq):
  return seq.translate(rev_table)[::-1]

def fetchSeq(qr,fasta,s_path):
  ''' s_path is the path of samtools  '''
 
#  print qstart,qend
  region = '%s:%d-%d'%(qr.chr, qr.start, qr.stop - 1)
  cmd_list = [s_path, "faidx", fasta, region]
  p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
  (std_out, std_err) = p.communicate()
  exit_code = p.returncode

  if exit_code == 0:
    seq = "".join(std_out.split('\n')[1:])
    if qr.strand == "-":
      seq = revcomp(seq)
    return seq.upper()
  else:
    raise Exception(205, std_err)

def fetchHistModSeq(qr, dbi):
  #[a,b)

  pointer=qr.start
  hseq=""

  interval_list=sorted(list(dbi.query(qr)))

  for line in interval_list:
    if not overlap(line,qr): 
      continue

    bstart=line.start
    bend=line.stop

    if bstart<pointer:
      hseq+="1"*(min(bend,qr.stop)-pointer)
    elif bstart>pointer:
      hseq+="0"*(bstart-pointer)+"1"*(min(bend,qr.stop)-bstart)
    else:
      hseq+="1"*(min(bend,qr.stop)-bstart)
    pointer=min(bend,qr.stop)

  if pointer!=qr.stop:
    hseq+="0"*(qr.stop-pointer)
  if qr.strand=="-":
    hseq=hseq[::-1]

  one_num = hseq.count("1")
  return hseq, one_num

def Generate_output_str(bed_pair):
  bed1=bed_pair[0]
  bed2=bed_pair[1]

  output_list=[]
  i = 0
  for s,sp,qrs in zip(args.species, [sp1,sp2], [bed1,bed2]):
    output_list.append("@"+s+"_"+"|".join(args.histone)+"_"+bed1.id)

    seq=fetchSeq(qrs,args.fasta_path.rstrip('/')+'/'+s+'.fa',args.s_path)
    if "N" in seq:
      output_list=[]
      return
    output_list.append(seq)
    if i == 0 :
      one_numbers = []
    for fsp_name in sp:
      output_list.append("+")
      hseq, one_num = fetchHistModSeq(qrs, fsp_name)
      if i == 0:
        one_numbers.append(one_num)
      if len(hseq)!=len(seq):
        # print >> sys.stderr, "In region pair " + bed1.id + " the sequence and epigenomic state strings have different lengths."
        raise Exception(203, "In region pair " + bed1.id + " the sequence and epigenomic state strings have different lengths.")
      output_list.append(str(hseq))

    i += 1
  one_name = "$$$" + "$".join([str(f) for f in one_numbers])
  output_list[0] += one_name
  output_list[4] += one_name

  if len(output_list)!=0:
    return "\n".join(output_list)


def Main():
  global args
  args=ParseArg()

  # Index peak files
  global sp1, sp2
  sp1, sp2=ReadHistones(args.bg)
 
  if len(sp1)!=len(sp2) or len(sp1)!=len(args.histone):
    print >> sys.stderr, "The number of histone marks must be identical for the two species. Provided histone mark names should match the peak files."
    sys.exit(204)

  fout=open(args.output,"w")

  # t0 = time()
  input_list=[]
  with open(args.q_region[0],"r") as fbed1, open(args.q_region[1],"r") as fbed2:
    while True:
      line1=fbed1.readline().strip().split()
      line2=fbed2.readline().strip().split()
      if len(line1) == 0:
        break
      if len(line2) == 0:
        break

      bed1=Bed(line1)
      bed2=Bed(line2)

      if bed1.chr != "chrX" and bed1.chr != "chrY" and (not bed1.chr.lstrip("chr").isdigit()):
        continue
      if bed2.chr != "chrX" and bed2.chr != "chrY" and (not bed2.chr.lstrip("chr").isdigit()):
        continue

      input_list.append((bed1, bed2))

  p = Pool(args.p_num)
  try:
    out_queue = p.map(Generate_output_str,input_list)
    p.close()
    p.join()
  except Exception as e:
    p.terminate()
    print >> sys.stderr, e.args
    sys.exit(e.args[0])

  # t1 = time()
  # print >>sys.stderr, "Time: " + str((t1 - t0) / 60)

  with open(args.output, "w") as fout:
    for item in out_queue:
      if item:
        print >>fout, item

Main()





  
