import sys, argparse
from xplib.Annotation import Bed
from xplib import DBI

def ParseArg():
  p=argparse.ArgumentParser(description="Prepare data for ribbon plot.")
  p.add_argument("input", type=str, help="Region bed file.")
  p.add_argument("-a", "--annotation", type=str, help="Annotation file.")
  p.add_argument("-e", "--ext_dis", type=int, help="Extension distance.")
  p.add_argument("-n", "--target_num", type=int, help="Target number.")
  p.add_argument("-o","--output",type=str, help="output file name.")
  # output format: regionId, geneId, distance. 
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args() 


def overlap(bed1,bed2):
  if (bed1.stop > bed2.start) and (bed1.start < bed2.stop):
    return True
  else:
    return False


def findNearbyGene(bed, dbi, ori_start, ori_stop, n):
  dis_list = []
  geneId_list = []

  interval_list = dbi.query(bed)

  for gene in interval_list:
    if not overlap(gene, bed): 
      continue
    dis = max(ori_start - gene.stop, gene.start - ori_stop)
    if dis < 0:
      dis = 0
    dis_list.append(dis)
    geneId_list.append(gene.id)
  
  gene_list_sort = sorted(zip(dis_list, geneId_list), reverse=True)[0:n]
  return gene_list_sort

def Main():
  args = ParseArg()

  anno = DBI.init(args.annotation, "bed")
  ext_dis = args.ext_dis
  target_num = args.target_num

  with open(args.input, "r") as fin, open(args.output, "w") as fout:
    for line in fin:
      bed_region = Bed(line.strip().split())
      mid_point = (bed_region.start + bed_region.stop) / 2
      ori_start = bed_region.start
      ori_stop = bed_region.stop
      bed_region.start = mid_point - ext_dis
      bed_region.stop = mid_point + ext_dis

      gene_list = findNearbyGene(bed_region, anno, ori_start, ori_stop, target_num)
      for gene in gene_list:
        print >> fout, "\t".join([bed_region.id, gene[1], str(gene[0])])
  
Main()




