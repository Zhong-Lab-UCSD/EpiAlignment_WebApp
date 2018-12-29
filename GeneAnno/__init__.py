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
  Cut the promoter regions and return a bed 6 list. All elements are strings.
  '''
  if trans_obj.strand == "+":
    start = trans_obj.start - cutUp
    stop = trans_obj.start + cutDown
  else:
    start = trans_obj.stop - cutDown
    stop = trans_obj.stop + cutUp
  return [trans_obj.chr, str(start), str(stop), trans_obj.ensid + "_" + trans_obj.transname, '0', trans_obj.strand]



