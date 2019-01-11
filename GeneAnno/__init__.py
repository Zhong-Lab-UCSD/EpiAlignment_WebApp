from copy import deepcopy

class Transcripts:
  def __init__(self, x):
    '''
    x: a list with 9 elements:
    chr, start, stop, genename, ensid, transname, genetype, strand, transtype.
    '''
    self.chr = x[0]
    self.start = int(x[1])
    self.stop = int(x[2])
    self.name = x[3]
    self.ensid = x[4]
    self.transname = x[5]
    self.genetype = x[6]
    self.strand = x[7]
    self.transtype = x[8]

def Construct_ensDict(fname):
  ensDict = {}
  with open(fname, "r") as fin:
    for line in fin:
      line = line.strip().split("\t")
      trans_obj = Transcripts(line)
      if trans_obj.genetype == "protein_coding" and trans_obj.transtype != "protein_coding":
        continue

      if trans_obj.ensid not in ensDict:
        ensDict[trans_obj.ensid] = []
      ensDict[trans_obj.ensid].append(trans_obj)
  return ensDict

def Construct_nameDict(fname):
  nameDict = {}
  with open(fname, "r") as fin:
    for line in fin:
      line = line.strip().split("\t")
      trans_obj = Transcripts(line)
      if trans_obj.genetype == "protein_coding" and trans_obj.transtype != "protein_coding":
        continue

      if trans_obj.name not in nameDict:
        nameDict[trans_obj.name] = []
      nameDict[trans_obj.name].append(trans_obj)
  return nameDict

def PromoterBed(trans_obj, cutUp = 0, cutDown = 0):
  '''
  Cut the promoter regions and return a bed 6 list.
  '''
  if trans_obj.strand == "+":
    start = trans_obj.start - cutUp
    stop = trans_obj.start + cutDown
  else:
    start = trans_obj.stop - cutDown
    stop = trans_obj.stop + cutUp
  return [trans_obj.chr, start, stop, trans_obj.ensid + "_" + trans_obj.name, '0', trans_obj.strand]

def Dist_overlap(region1, region2, N):
  start1 = region1[1]
  stop1 = region1[2]
  start2 = region2[1]
  stop2 = region2[2]

  if min((stop1 - start2), (stop2 - start1)) >= N:
    return True
  else:
    return False


def PromoterMerge(gene, trans_dict, cutUp = 0, cutDown = 0):
  '''
  Merge overlapping promoters of the same gene.
  '''
  N = float(cutUp + cutDown) * 0.8
  res_list = []
  if gene not in trans_dict:
    return res_list
  promoter_list = [PromoterBed(trans_obj, cutUp, cutDown) for trans_obj in trans_dict[gene]]
  # sort tmp_list by chromosome names and starts.
  promoter_list_sort = sorted(promoter_list, key=lambda x: x[1])
  for p in promoter_list_sort:
    if len(res_list) == 0:
      res_list.append(deepcopy(p))
    if Dist_overlap(p, res_list[-1], N):
      # Update start
      res_list[-1][1] = min(p[1], res_list[-1][1])
      # Update stop
      res_list[-1][2] = max(p[2], res_list[-1][2])
    else:
      # change region name
      res_list[-1][1] = str(res_list[-1][1])
      res_list[-1][2] = str(res_list[-1][2])
      res_list[-1][3] += res_list[-1][1][-4:]
      res_list.append(deepcopy(p))
  # The last one
  res_list[-1][1] = str(res_list[-1][1])
  res_list[-1][2] = str(res_list[-1][2])
  res_list[-1][3] += res_list[-1][1][-4:]
  return res_list



