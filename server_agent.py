import subprocess
from subprocess import Popen, PIPE
from numpy import percentile, mean
from scipy.stats import norm
from itertools import izip
from GeneAnno import *
import json
import shutil
import sys
import os
import re

ensembl_regexp = 'ENS[A-Z]+[0-9]{11}'
WARNING_SIZE = 10000

def ParseJson():
  '''
  Parse the json string passed by nodejs from stdin.
  return: a dictionary with three keys: body, files, runid (output folder name).
  '''
  lines = sys.stdin.readlines()
  json_dict = json.loads(lines[0])
  return json_dict


def MoveUploadFiles(dest_name, fdict):
  '''
  Move uploaded files to the output folder.
  dest_name: output folder name.
  fdict: the value of "files" of the json dictionary.
  '''
  for key in fdict:
    for ufile in fdict[key]:
      shutil.move(ufile["path"], dest_name + ufile["filename"])


def CheckFileLength(file1, file2):
  with open(file1, "r") as fin1, open(file2, "r") as fin2:
    if len(fin1.readlines()) == len(fin2.readlines()):
      return True
    else:
      print >> sys.stderr, "[EpiAlignment]The two input files have different numbers of regions."
      sys.exit(202)


def CheckLineType(line, oldType = None):
  line_len = len(line)
  if (line_len == 6 and oldType != "name"):
    # the line has 6 fields.
    return "bed"
  elif (line_len == 1 and oldType != "bed"):
    # the line has 1 fields
    return "name"
  else:
    if line_len == 1 or line_len == 6:
      print >> sys.stderr, "[EpiAlignment]The format of your input file is not consistent. Please use all gene names or all BED6 format."
    else:
      print >> sys.stderr, "[EpiAlignment]Input files have to be bed6 files (6 columns) or genelists (1 columns)."
    sys.exit(201)

def StripDigits(qstr):
  '''
  Remove all digits in qstr.
  '''
  return ''.join([i for i in qstr if not i.isdigit()])


def ParsePeaks(of_name, json_dict, runid):
  '''
  Creak peak files.
  preset_data: a list with a pair of peak file id.
  '''
  json_files = json_dict["files"]
  db_fname = "html/assets/experimentDict.json"
  if "encodeData" in json_dict["body"]:
    # preset data select
    preset_data = json_dict["body"]["encodeData"]
    with open(db_fname, "r") as fdb:
      data_json = json.load(fdb)
    peak1 = data_json[preset_data[0]]['peak_file']
    peak2 = data_json[preset_data[1]]['peak_file']
  else:
    # upload peak files
    if "speciesPeak1[]" not in json_files or "speciesPeak2[]" not in json_files :
      print >> sys.stderr, "[EpiAlignment]No peak files found!"
      sys.exit(207)
    peak1 = of_name + json_files["speciesPeak1[]"][0]["filename"]
    peak2 = of_name + json_files["speciesPeak2[]"][0]["filename"]

  with open(of_name + "peaks_" + runid, "w") as fpeak:
    print >> fpeak, "@species1"
    print >> fpeak, peak1
    print >> fpeak, "@species2"
    print >> fpeak, peak2


def parseBed(bedFields, nameDict, warningSize=None, nameIndex=3, strandIndex=5):
  try:
    suffix = 1
    if bedFields[nameIndex] in nameDict:
      suffix = nameDict[bedFields[nameIndex]]
      while (bedFields[nameIndex] + '_' + str(suffix)) in nameDict:
        suffix += 1
      bedFields[nameIndex] = bedFields[nameIndex] + '_' + str(suffix)
    nameDict[bedFields[nameIndex]] = suffix + 1
    if bedFields[strandIndex] == "+" or bedFields[strandIndex] == "." or bedFields[strandIndex] == "1":
      bedFields[strandIndex] = "+"
    elif bedFields[strandIndex] == "-" or bedFields[strandIndex] == "0":
      bedFields[strandIndex] = "-"
    else:
      raise Exception('strand is not "+", "-" or "." ("." will be considered as "+").')
    if warningSize is not None and int(bedFields[2]) - int(bedFields[1]) > int(warningSize):
      print >> sys.stderr, ('[EpiAlignment]Long query region detected (>' +
                            str(warningSize) + 'bp). ' +
                            'EpiAlignment is designed to match ' +
                            'medium-sized functional genomic elements to '+
                            'a best hit within the long target region. ' +
                            'The biological insight from results of such ' +
                            'long query regions may be limited. ' +
                            '(Please see the manual for details.) ' +
                            'In addition, the run time may be very long.')
    return bedFields
  except (IndexError, TypeError):
    raise Exception('Not a BED6 format.')
      

def FileOrTextarea(textarea_input, json_files, key, of_name, runid, warningSize = None):
  '''
  Determine if input was pasted into the textarea or uploaded as a file.
  textarea_input: a string. Text in textarea.
  json_files: json "files" value.
  key: key value of the file.
  align_mode: search mode, promoter or enhancer.
  of_name: output folder name.
  return: a string and file type. If a file was uploaded, simply return file name. If data was pasted into the textarea,
  write the data into a new file.
  '''
  lineType = None
  regionNameDict = dict()
  if key in json_files:
    # uploaded file
    fname = of_name + json_files[key][0]["filename"]
    fOutName = of_name + key + "_" + runid
    with open(fname, 'r') as fIn, open(fOutName, 'w') as fOut:
      lineNum = 0
      for line in fIn:
        try:
          lineNum += 1
          line = line.strip().split()
          lineType = CheckLineType(line, lineType)
          if lineType == 'bed':
            parseBed(line, regionNameDict, warningSize)
          print >> fOut, '\t'.join(line)
        except Exception as err:
          print >> sys.stderr, '[EpiAlignment]Skipping line #' + str(
              lineNum) + ': ' + err.message
  elif textarea_input != "":
    # paste data
    fOutName = of_name + key + "_" + runid
    with open(fOutName, "w") as fOut:
      lineNum = 0
      lines = textarea_input.rstrip("\n").split('\n')
      for line in lines:
        try:
          lineNum += 1
          line = line.strip().split()
          lineType = CheckLineType(line, lineType)
          if lineType == 'bed':
            parseBed(line, regionNameDict, warningSize)
          print >> fOut, '\t'.join(line)
        except Exception as err:
          print >> sys.stderr, '[EpiAlignment]Skipping line #' + str(
              lineNum) + ': ' + err.message
  else:
    # no data provided.
    return "", ""

  return fOutName, lineType

########################
## Mode 1: genelist   ##
########################
def Cons_transDict(gene_name, sp_name):
  if re.match(ensembl_regexp, gene_name):
    transDict = Construct_ensDict("Annotation/AnnotationFiles/" + sp_name + "_transcript.clean")
  else:
    transDict = Construct_nameDict("Annotation/AnnotationFiles/" + sp_name + "_transcript.clean")
  return transDict


def Cons_transList(input1, intype1, promoterUp, promoterDown, sp):
  trans_list1 = []
  with open(input1, "r") as fin1:
    if intype1 == "bed":
      trans_list1 = [line.strip().split() for line in fin1]
    elif intype1 == "name":
      i = 0
      for line in fin1:
        line = line.strip()
        if i == 0:
          transDict1 = Cons_transDict(line, sp)
          i += 1
        if line in transDict1:
          trans_list1 += PromoterMerge(line, transDict1, promoterUp, promoterDown)
        else:
          print >> sys.stderr, "[EpiAlignment]The gene " + line + " was not found in " + sp
  return trans_list1


def PairCutPromoter(input1, input2, intype1, intype2, promoterUp, promoterDown, genAssem):
  trans_list1 = Cons_transList(input1, intype1, promoterUp, promoterDown, genAssem[0])
  trans_list2 = Cons_transList(input2, intype2, promoterUp, promoterDown, genAssem[1])

  with open(input1 + ".bed", "w") as fout1, open(input2 + ".bed", "w") as fout2:
    i = 0
    for region1 in trans_list1:
      for region2 in trans_list2:
        region_name = region1[3] + "[===]" + region2[3]
        print >> fout1, "\t".join(region1[0:3] + [region_name] + region1[4:])
        print >> fout2, "\t".join(region2[0:3] + [region_name] + region2[4:])
        i += 1

  if i > 10000:
    print >> sys.stderr, "[EpiAlignment]Too many regions..."
    sys.exit(210)

  return input1 + ".bed", input2 + ".bed"

def PairCutEnhancer(input1, input2, promoterUp, promoterDown, genAssem):
  '''
  Pair promoters (multiple) and bed regions in the enhancer mode when
  query regions are provided by gene names and target regions are provided as bed regions.
  '''
  i = 0
  with open(input1, "r") as fin1, open(input2, "r") as fin2, \
    open(input1 + ".bed", "w") as fout1, open(input2 + ".bed", "w") as fout2:
    for name, bed in izip(fin1, fin2):
      name = name.strip()
      bed = bed.strip().split()
      if i == 0:
        transDict = Cons_transDict(name, genAssem[0])
        i += 1
      if name in transDict:
        trans_list = PromoterMerge(name, transDict, promoterUp, promoterDown)
      else:
        print >> sys.stderr, "[EpiAlignment]Gene %s is not found and skipped."%(name)
      for region in trans_list:
        region_name = region[3] + "Vs" + bed[3]
        print >> fout1, "\t".join(region[0:3] + [region_name] + region[4:])
        print >> fout2, "\t".join(bed[0:3] + [region_name] + bed[4:])
  return input1 + ".bed", input2 + ".bed"
      
    
#######################
## Mode 3: cluster   ##
#######################
def GenesInCluster(cluster_id, sp, of_name):
  '''
  Extract all genes in the selected cluster.
  cluster_id: the ID of the cluster. 
  sp: genome assembly name.
  return: cluster_genes, a list of gene ensembl ids.
  '''
  cluster_genes = set()
  cfname =  "Annotation/AnnotationFiles/" + sp + "_clusters"
  with open(cfname, "r") as cfin:
    for line in cfin:
      line = line.strip().split()
      if line[2] == cluster_id:
        cluster_genes.add(line[0])
  return list(cluster_genes)

def PairCutCluster(input1, intype1, cluster_id, promoterUp, promoterDown, genAssem, runid, of_name):
  '''
  Find genes in the selected cluster.
  '''
  # Extract genes in the cluster
  cluster_genes2 = GenesInCluster(cluster_id, genAssem[1], of_name)
  transDict2 = Cons_transDict(cluster_genes2[0], genAssem[1])
  trans_list2 = []
  for gene in cluster_genes2:
    if gene in transDict2:
      trans_list2 += PromoterMerge(gene, transDict2, promoterUp, promoterDown)
  
  fname2 = of_name + cluster_id + runid + "_2.bed"
  if input1 != "":
    # uploaded file
    trans_list1 = Cons_transList(input1, intype1, promoterUp, promoterDown, genAssem[0])
    with open(input1 + ".bed", "w") as fout1, open(fname2, "w") as fout2:
      for region1 in trans_list1:
        for region2 in trans_list2:
          region_name = region1[3] + "[===]" + region2[3]
          print >> fout1, "\t".join(region1[0:3] + [region_name] + region1[4:])
          print >> fout2, "\t".join(region2[0:3] + [region_name] + region2[4:])
    return input1 + ".bed", fname2

  else:
    # Input1 is empty
    fname1 = of_name + cluster_id + runid + "_1.bed"
    cluster_genes1 = GenesInCluster(cluster_id, genAssem[0], of_name)
    transDict1 = Cons_transDict(cluster_genes1[0], genAssem[0])
    trans_list1 = []
    for gene in cluster_genes1:
      if gene in transDict1:
        trans_list1 += PromoterMerge(gene, transDict1, promoterUp, promoterDown)

    with open(fname1, "w") as fout1, open(fname2, "w") as fout2:
      for region1 in trans_list1:
        for region2 in trans_list2:
          region_name = region1[3] + "[===]" + region2[3]
          print >> fout1, "\t".join(region1[0:3] + [region_name] + region1[4:])
          print >> fout2, "\t".join(region2[0:3] + [region_name] + region2[4:])

    return fname1, fname2 



#######################
## Mode 4: liftOver  ##
#######################
def GeneNameToBed(input1, promoterUp, promoterDown, genAssem):
  '''
  This function convert gene names to promoter beds tn the one-vs-one mode
  when the query regions are provided by gene names.
  '''
  i = 0
  with open(input1, "r") as fin1, open(input1 + ".bed", "w") as fout1:
    for name in fin1:
      name = name.strip()
      if i == 0:
        transDict = Cons_transDict(name, genAssem[0])
        i += 1
      if name in transDict:
        trans_list = PromoterMerge(name, transDict, promoterUp, promoterDown)
      else:
        print >> sys.stderr, "[EpiAlignment]Gene %s is not found and skipped."%(name)
      for region in trans_list:
        print >> fout1, "\t".join(region)
  return input1 + ".bed"


def ExtendBed(fname, enhUp, enhDown):
  '''
  Extend the input bed file for liftOver.
  return: extbed, name of the file with extended regions.
  '''
  with open(fname, "r") as fin, open(fname + ".extend", "w") as fout:
    for line in fin:
      line = line.strip().split()
      if line[5] == "+":
        line[1] = str(int(line[1]) - enhUp)
        line[2] = str(int(line[2]) + enhDown)
      else:
        line[1] = str(int(line[1]) - enhDown)
        line[2] = str(int(line[2]) + enhUp)
      print >> fout, "\t".join(line)

  return fname + ".extend"


def LiftOver(input_bed, genAssem):
  '''
  Call liftOver to remap regions.
  '''
  chain_name = genAssem[0] + "To" + genAssem[1].capitalize() + ".over.chain"
  cmd_list = ["liftOver", input_bed, "Annotation/AnnotationFiles/" + chain_name, input_bed + ".lift", input_bed + ".unlift" ,"-minMatch=0.1"]
  p = subprocess.Popen(cmd_list, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  (std_out, std_err) = p.communicate()
  exit_code = p.returncode

  if exit_code != 0:
    print >> sys.stderr, "[EpiAlignment]Failed to generate the input file. liftOver exited with code: " + str(exit_code)
    sys.exit(exit_code)
  return input_bed + ".lift"

def RemoveNonlift(input_bed, lift_bed):
  '''
  Remove non-remappable lines from the input_bed.
  '''
  i = 0
  with open(input_bed, "r") as fin1, open(lift_bed, "r") as fin2, open(input_bed + ".clean", "w") as fout:
    while True:
      line2 = fin2.readline().strip()
      if line2 == "":
        if i == 0:
          # lift_bed is empty
          print >> sys.stderr, "[EpiAlignment]None of the input regions were remappable."
          sys.exit(209)
        break
      i += 1
      line2 = line2.split()

      while True:
        line1 = fin1.readline().strip().split()
        if line1[3] == line2[3]:
          print >> fout, "\t".join(line1) 
          break
        if len(line1) == 0:
          break
  return input_bed + ".clean"

##################
## Create beds  ##
##################
def CheckNumber(key, json_body):
  if key in json_body:
    return int(json_body[key])
  else:
    return 0

def CreateInputBeds(of_name, json_dict, runid):
  '''
  Create a pair of bed files.
  of_name: output folder name.
  json_dict: the json dictionary.
  return: names of the two bed files.
  '''
  # Common variables
  searchMode = json_dict["body"]["searchRegionMode"]
  alignMode = json_dict["body"]["alignMode"]
  genAssem = json_dict["body"]["genomeAssembly"]
  promoterUp = CheckNumber("promoterUp", json_dict["body"])
  promoterDown = CheckNumber("promoterDown", json_dict["body"])
  enhancerUp = CheckNumber("enhancerUp", json_dict["body"])
  enhancerDown = CheckNumber("enhancerDown", json_dict["body"])

  # Is input1 a file or a pasted text?
  input1, intype1 = FileOrTextarea(json_dict["body"]["speciesText"][0], json_dict["files"], "speciesInput1", of_name, runid, WARNING_SIZE)
  if input1 == "" and (not json_dict["body"]["searchRegionMode"] == "genecluster"):
    print >> sys.stderr, "[EpiAlignment]No input regions provided."
    sys.exit(200)

  if searchMode == "genomeregion":
    # Mode 1: define search regions with bed files or gene lists.
    # Is input2 a file or a pasted text?
    input2, intype2 = FileOrTextarea(json_dict["body"]["speciesText"][1], json_dict["files"], "speciesInput2", of_name, runid)
    if alignMode == "enhancer":
      if CheckFileLength(input1, input2):
        if intype2 == "name":
          print >> sys.stderr, "[EpiAlignment]In one-vs-one mode, only query regions can be defined with gene names. \
            Target regions need to be provided in BED6 format."
          sys.exit(201)
        else:
          if intype1 == "bed":
            return input1, input2, intype1, intype2
          elif intype1 == "name":
            bed1, bed2 = PairCutEnhancer(input1, input2, promoterUp, promoterDown, genAssem)
    else:
      bed1, bed2 = PairCutPromoter(input1, input2, intype1, intype2, promoterUp, promoterDown, genAssem)
    return bed1, bed2, intype1, intype2

  else:    
    if searchMode == "genetype" and alignMode == "promoter":
      # Mode 2: search the promoter regions of a specific type of gene.
      # Mode removed temporarily.
      pass

    elif searchMode == "genecluster" and alignMode == "promoter":
      # Mode 3: search a specific gene cluster.
      cluster_id = json_dict["body"]["clusters"]
      bed1, bed2 = PairCutCluster(input1, intype1, cluster_id, promoterUp, promoterDown, genAssem, runid, of_name)

    elif searchMode == "homoregion" and alignMode == "enhancer":
      # Mode 4 (enhancer mode 2): use homologous regions. 
      # species must be different!
      ori_input = input1
      if intype1 == "name":
        ori_input = GeneNameToBed(input1, promoterUp, promoterDown, genAssem)
      if StripDigits(genAssem[0]) == StripDigits(genAssem[1]):
        print >> sys.stderr, "[EpiAlignment]The two species must be different to use this mode."
        sys.exit(208)
      # Extend the input bed file1. extbed: extended bed file name.
      extbed = ExtendBed(ori_input, enhancerUp, enhancerDown)
      # LiftOver
      liftbed = LiftOver(extbed, genAssem)
      # Remove non-remappable regions. Return a pair of bed file names.
      cleanbed = RemoveNonlift(ori_input, liftbed)
      #os.remove(extbed)
      #os.remove(extbed + ".unlift")
      bed1 = cleanbed
      bed2 = liftbed

    return bed1, bed2, intype1, ""


def BedToFa(bed1, bed2, out_folder, sp_list, runid):
  '''
  Run InputToFastq_bed2.py to convert the bed file pair to a fastq-like file.
  '''
  # epiName_list = epiName.split(",")
  cmd_list = ["python", "InputToFastq_bed2.py", bed1, bed2, "-s"] + sp_list +\
  ["--bg", out_folder + "peaks_" + runid] +\
  ["--histone", "epi"] +\
  ["-p", "20"] +\
  ["-o", out_folder + "Input_" + runid]

  p = Popen(cmd_list, stderr=PIPE)
  (std_out, std_err) = p.communicate()
  exit_code = p.returncode

  if exit_code != 0:
    print >> sys.stderr, "[EpiAlignment]" + std_err + " Exit code: " + str(exit_code)
    sys.exit(exit_code)
  # Check if the input file is empty.
  Input_fa_name = out_folder + "Input_" + runid
  if os.stat(Input_fa_name).st_size == 0:
    print >> sys.stderr, "[EpiAlignment]Failed to generate the input file for EpiAlignment. Please check whether your genomic regions/gene names match the genome assemblies."
    sys.exit(214)


def InputParas(of_name, json_body, runid):
  '''
  Create input parameter file for EpiAlignment.
  '''
  # Check if parameters are not positive.
  if "seqweight" in json_body:
    seqweight = json_body["seqweight"]
  else:
    seqweight = "1"
  parak_list = [float(x) for x in json_body["parak"].split(",")]
  if float(json_body["paras"]) <= 0 or float(json_body["paramu"]) <= 0 or min(parak_list) <= 0:
    print >> sys.stderr, "[EpiAlignment]Parameters should be positive values."
    sys.exit(206)

  seq_pi_list = [float(json_body["piA"]), float(json_body["piC"]), float(json_body["piG"]), float(json_body["piT"])]
  pi_list1 = [float(k) for k in json_body["pi1"].split(",")]
  weight_list = [float(w) for w in json_body["epiweight"].split(",")]
  para_list = seq_pi_list + pi_list1
  if min(para_list) <= 0 or max(para_list) >= 1:
    print >> sys.stderr, "[EpiAlignment]Equilibrium probabilities (pi) must be values between 0 and 1."
    sys.exit(206)

  # Create parameter file.
  with open(of_name + "parameters_" + runid, "w") as fpara:
    print >> fpara, json_body["paras"]
    print >> fpara, json_body["paramu"]
    print >> fpara, "\n".join(json_body["parak"].split(","))
    print >> fpara, "A:" + json_body["piA"] + "\t" + "C:" + json_body["piC"] + "\t" +\
    "G:" + json_body["piG"] + "\t" + "T:" + json_body["piT"]
    pi_list0 = [1 - k for k in pi_list1]
    for p0, p1 in zip(pi_list0, pi_list1):
      print >> fpara, "0:" + str(p0) + "\t" + "1:" + str(p1)
    weights = "\t".join(json_body["epiweight"].split(","))
    print >> fpara, seqweight + "\t" + weights

  # If epi weight is not 0, create another parameter file for sequence-only alignment.
  if json_body["epiweight"] != "0":
    with open(of_name + "parameters_seq_" + runid, "w") as fseq_para:
      print >> fseq_para, json_body["paras"]
      print >> fseq_para, json_body["paramu"]
      print >> fseq_para, "\n".join(json_body["parak"].split(","))
      print >> fseq_para, "A:" + json_body["piA"] + "\t" + "C:" + json_body["piC"] + "\t" +\
      "G:" + json_body["piG"] + "\t" + "T:" + json_body["piT"]
      pi_list0 = [1 - k for k in pi_list1]
      for p0, p1 in zip(pi_list0, pi_list1):
        print >> fseq_para, "0:" + str(p0) + "\t" + "1:" + str(p1)
      print >> fseq_para, "1" + "\t" + "0"


def ExeEpiAlignment(alignMode, searchRegionMode, bed1, bed2, genAssem, of_name, runid):
  '''
  Execute EpiAlignment
  '''
  seq_stat = os.path.isfile(of_name + "parameters_seq_" + runid)
  cmd_list = ["python3", "EpiAlignment_3.py", of_name + "Input_" + runid] +\
    ["-e", of_name + "parameters_" + runid] +\
    ["-p", "140"] +\
    ["-o", of_name + "epialign_res_" + runid]

  cmd_list_seq = ["python3", "EpiAlignment_3.py", of_name + "Input_" + runid] +\
    ["-e", of_name + "parameters_seq_" + runid] +\
    ["-p", "140"] +\
    ["-o", of_name + "seqalign_res_" + runid]

  # Fetch gene Ids.
  cmd_list_gene1 = ["python", "EnhancerOverlappingGenes.py", bed1] + \
    ["-a", "Annotation/AnnotationFiles/genes/" + genAssem[0] + ".genes.ensembl.sorted.txt" ] + \
    ["-e", "50000"] + \
    ["-n", "5"] + \
    ["-o", of_name + "QueryRNA_" + runid]
  
  cmd_list_gene2 = ["python", "EnhancerOverlappingGenes.py", bed2] + \
    ["-a", "Annotation/AnnotationFiles/genes/" + genAssem[1] + ".genes.ensembl.sorted.txt" ] + \
    ["-e", "50000"] + \
    ["-n", "5"] + \
    ["-o", of_name + "TargetRNA_" + runid]

  if alignMode == "promoter":
    p_epi = Popen(cmd_list, stderr=PIPE)
    if seq_stat:
      p_seq = Popen(cmd_list_seq, stderr=PIPE)
    (std_out_epi, std_err_epi) = p_epi.communicate()
    exit_code_epi = p_epi.returncode
    if exit_code_epi != 0:
      print >> sys.stderr, "[EpiAlignment]Failed to align regions. Exit code: " + str(exit_code_epi)
      sys.exit(exit_code_epi)

    if seq_stat:
      (std_out_seq, std_err_seq) = p_seq.communicate()
      exit_code_seq = p_seq.returncode
      if exit_code_seq != 0:
        print >> sys.stderr, "[EpiAlignment]Failed to align regions. Exit code: " + str(exit_code_seq)
        sys.exit(exit_code_seq)

  elif alignMode == "enhancer":
    cmd_list += ["-O", of_name + "epi_scores_" + runid]
    cmd_list_seq += ["-O", of_name + "seq_scores_" + runid]

    p_epi = Popen(cmd_list, stderr=PIPE)
    if seq_stat:
      p_seq = Popen(cmd_list_seq, stderr=PIPE)
    # Fetch overlapping genes.
    p_gene1 = Popen(cmd_list_gene1, stderr=PIPE)
    p_gene2 = Popen(cmd_list_gene2, stderr=PIPE)

    (std_out_gene1, std_err_gene1) = p_gene1.communicate()
    (std_out_gene2, std_err_gene2) = p_gene2.communicate()
    exit_code_gene1 = p_gene1.returncode
    exit_code_gene2 = p_gene2.returncode

    if exit_code_gene1 != 0 or exit_code_gene2 != 0:
      print >> sys.stderr, "[EpiAlignment]Error occurred when fetching expression data."

    (std_out_epi, std_err_epi) = p_epi.communicate()
    exit_code_epi = p_epi.returncode
    if exit_code_epi != 0:
      print >> sys.stderr, "[EpiAlignment]Failed to align regions. Exit code: " + str(exit_code_epi)
      sys.exit(exit_code_epi)

    if seq_stat:
      (std_out_seq, std_err_seq) = p_seq.communicate()
      exit_code_seq = p_seq.returncode
      if exit_code_seq != 0:
        print >> sys.stderr, "[EpiAlignment]Failed to align regions. Exit code: " + str(exit_code_seq)
        sys.exit(exit_code_seq)

###################
## Parse results ##
###################

def BedDict(fname):
  bed_dict = {}
  with open(fname, "r") as fin:
    for line in fin:
      line = line.strip().split()
      bed_dict[line[3]] = line[0:3] + [line[5]]
  return bed_dict

def TargetRegion(bed_list, hit_start, hit_stop):
  if bed_list[3] == "+":
    start = int(bed_list[1]) + int(hit_start)
    stop = int(bed_list[1]) + int(hit_stop)
  elif bed_list[3] == "-":
    start = int(bed_list[2]) - int(hit_stop)
    stop = int(bed_list[2]) - int(hit_start)
  return bed_list[0] + ":" + str(start) + "-" + str(stop) + "(" + bed_list[3] + ")"

def ConcateBed(coor_list):
  return coor_list[0] + ":" + str(coor_list[1]) + "-" + str(coor_list[2]) + "(" + coor_list[3] + ")"

def InitJsonObj(ind, pair_name, bed_dict1, bed_dict2, line_epi, line_seq, one_num):
  '''
  Initialize a json object: index; all locations including
  region1, region2, queryLength,
  targetE, target S; epi-hit scoreE and seq-hit scoreS.
  ifseq: If sequence-only alignment has been done. 
  '''
  query_len = int(bed_dict1[pair_name][2]) - int(bed_dict1[pair_name][1])
  json_obj = {"index":ind, "region1": ConcateBed(bed_dict1[pair_name]), "region2": ConcateBed(bed_dict2[pair_name]),\
    "queryLength":query_len, \
    "scoreE": float(line_epi[1]) * 1000 / query_len, "targetE": TargetRegion(bed_dict2[pair_name], line_epi[5], line_epi[6]), \
    "scoreS": ".", "targetS": ".", "shifted": ".", "oneNum": one_num}

  if line_seq:
    json_obj["scoreS"] = float(line_seq[1]) * 1000 / query_len
    json_obj["targetS"] = TargetRegion(bed_dict2[pair_name], line_seq[5], line_seq[6])
    # shifted or not
    if abs(int(line_epi[6]) - int(line_seq[6])) > json_obj["queryLength"]:
      json_obj["shifted"] = "Y"
    else:
      json_obj["shifted"] = "N"

  return json_obj


def RegionName(json_obj, pair_name, intype1, intype2, alignMode):
  name_dict = {"region_name1": ".", "region_name2": ".", "ensID1": ".", "ensID2": ".", "transID1":".", "transID2": "."}
  json_obj.update(name_dict)
  if (alignMode == "enhancer"):
    # a pair of bed.
    json_obj["region_name1"] = pair_name
  else:
    pair_name = pair_name.split("[===]")
    if intype1 == "bed":
      json_obj["region_name1"] = pair_name[0]
    else:
      json_obj["ensID1"] = pair_name[0].split("_")[0]
      json_obj["transID1"] = pair_name[0].split("_")[1]
    if intype2 == "bed":
      json_obj["region_name2"] = pair_name[1]
    else:
      json_obj["ensID2"] = pair_name[1].split("_")[0]
      json_obj["transID2"] = pair_name[1].split("_")[1]

def ExtractScore(scores, pos, ext_dis):
  start = max(0, pos - ext_dis)
  stop = pos + ext_dis
  return max(scores[start:stop])

def Signal_to_Noise(scores, bin, query_len):
  '''
  Fragmentize the score list into bins with length bin. 
  Compute max and min in each bin. 
  Return medians of maximum and minimum. 
  '''
  i = query_len
  max_values = []
  min_values = []
  while i + bin < len(scores):
    max_values.append(max(scores[i:i + bin]))
    min_values.append(min(scores[i:i + bin]))
    i = i + bin
  max_values.append(max(scores[i:]))
  min_values.append(min(scores[i:]))
  return percentile(max_values, 75), percentile(min_values, 25)


def snCalculater(signal, mid_point, half_noise):
  if signal:
    return (signal - mid_point) / half_noise
  return "."

def SeqBg(s, mu, alignMode):
  seq_dict = {"backgroundMean": ".", "backgroundSd":".", "backgroundMedian": ".", "backgroundQ75": ".", \
    "backgroundQ25":".", "orthoMedian": ".", "orthoQ75": ".", "orthoQ25": "."}
  s = str(round(float(s), 2))
  mu = str(round(float(mu), 2))
  if alignMode == "enhancer":
    bg_anno = "Annotation/AnnotationFiles/enhancerBackground.txt"
  else:
    bg_anno = "Annotation/AnnotationFiles/promoterBackground.txt"
  with open(bg_anno, "r") as fin:
    for line in fin:
      line = line.strip().split()
      if line[0] == s and line[1] == mu:
        seq_dict["backgroundMean"] = float(line[2])
        seq_dict["backgroundSd"] = float(line[3])
        seq_dict["backgroundMedian"] = float(line[4])
        seq_dict["backgroundQ25"] = float(line[5])
        seq_dict["backgroundQ75"] = float(line[6])
        seq_dict["orthoMedian"] = float(line[7])
        seq_dict["orthoQ25"] = float(line[8])
        seq_dict["orthoQ75"] = float(line[9])
        break
  return seq_dict

def FitNorm(signal, mean_value, sd_value):
  if signal:
    return 1 - norm(mean_value, sd_value).cdf(signal)
  return "."


def SequenceEvaluation(json_obj, line_epi, line_seq, epiScore, seqScore, s, mu, seq_bg):
  '''
  In json object, "shifted" has three possible values: Y, N, .
  The last value means that only sequence-only alignment was performed.
  '''
  query_len = json_obj["queryLength"]
  norm_factor = 1000.0 / query_len
  seqEval_dict = {"scoreS2": ".", "scoreE2": "."}
  # Extract scores.
  s1 = json_obj["scoreS"]
  s2 = None
  if json_obj["shifted"] == "Y":
    # extract epi-score for the seq-hit
    e2 = ExtractScore(epiScore, int(line_seq[6]), 50) * norm_factor
    # extract seq-score for the epi-hit
    s2 = ExtractScore(seqScore, int(line_epi[6]), 50) * norm_factor
    seqEval_dict["scoreE2"] = e2
    seqEval_dict["scoreS2"] = s2
  elif json_obj["shifted"] == "N":
    seqEval_dict["scoreE2"] = json_obj["scoreE"]
    seqEval_dict["scoreS2"] = json_obj["scoreS"]
    s2 = s1
  elif json_obj["shifted"] == ".":
    s2 = json_obj["scoreE"]
    s1 = None
  # Orthologous and Background
  #seqEval_dict["bgPvalueS"] = FitNorm(s1, seq_bg["backgroundMean"], seq_bg["backgroundSd"])
  #seqEval_dict["bgPvalueE"] = FitNorm(s2, seq_bg["backgroundMean"], seq_bg["backgroundSd"])
  # SignalToNoise ratio
  if json_obj["shifted"] != ".":
    upper, lower = Signal_to_Noise(seqScore, 500, query_len)
  else:
    upper, lower = Signal_to_Noise(epiScore, 500, query_len)
  upper = upper * norm_factor
  lower = lower * norm_factor
  seqEval_dict["signalToNoise"] = {"upperBound": upper, "lowerBound": lower} 
  half_noise = (upper - lower) / 2
  mid_point = (upper + lower) / 2
  seqEval_dict["signalToNoise"]["snS"] = snCalculater(s1, mid_point, half_noise)
  seqEval_dict["signalToNoise"]["snE"] = snCalculater(s2, mid_point, half_noise)
  
  json_obj.update(seqEval_dict)

def WriteFinalResult(json_obj, fout, alignMode):
  if alignMode == "enhancer":
    if json_obj["index"] == 1:
      print >> fout, "\t".join(["Index", "Query_region_name", "Query_gene", "Query_transcript", "Query_coordinate",\
        "Target_region_name", "Target_gene", "Target_transcript", "Target_coordinate",\
        "EpiAlign_target", "EpiAlignHit_epiScore", "EpiAlignHit_seqScore", "EpiAlignHit_SNR", \
        "SeqOnly_target", "SeqOnlyHit_epiScore", "SeqOnlyHit_seqScore", "SeqOnlyHit_SNR", "HitShifted"])

    print >> fout, "\t".join([str(f) for f in [json_obj["index"],\
      json_obj["region_name1"], json_obj["ensID1"], json_obj["transID1"], json_obj["region1"],\
      json_obj["region_name2"], json_obj["ensID2"], json_obj["transID2"], json_obj["region2"],\
      json_obj["targetE"], json_obj["scoreE"], json_obj["scoreS2"], json_obj["signalToNoise"]["snE"], \
      json_obj["targetS"], json_obj["scoreE2"], json_obj["scoreS"], json_obj["signalToNoise"]["snS"], json_obj["shifted"] ] ])

  elif alignMode == "promoter":
    if json_obj["index"] == 1:
      print >> fout, "\t".join(["Index", "Query_region_name", "Query_gene", "Query_transcript", "Query_coordinate",\
        "Target_region_name", "Target_gene", "Target_transcript", "Target_coordinate",\
        "EpiAlign_target", "epiScore",  \
        "SeqOnly_target", "seqScore"])

    print >> fout, "\t".join([str(f) for f in [json_obj["index"],\
      json_obj["region_name1"], json_obj["ensID1"], json_obj["transID1"], json_obj["region1"],\
      json_obj["region_name2"], json_obj["ensID2"], json_obj["transID2"], json_obj["region2"],\
      json_obj["targetE"], json_obj["scoreE"], \
      json_obj["targetS"], json_obj["scoreS"] ] ])


def ParseAlignResults(bed1, bed2, intype1, intype2, alignMode, searchRegionMode, of_name, runid, s, mu):
  '''
  Parse alignment results.
  bed1, bed2: the two bed files used for generating input file.
  intype: "bed", "name" or an empty string.
  of_name: output folder.
  runid: runid.
  return: None. This function will write a json object to a file.
  json object items: index, region1(chr:start:stop), region2,
  queryLength, scoreE, targetE, scoreEalt, scoreS, targetS, scoreSalt,
  sequenceEval(signalToNoise, background, orthologous)
  region_name1, region_name2, ensID1, ensID2, transID1, transID2.
  '''
  bed_dict1 = BedDict(bed1)
  bed_dict2 = BedDict(bed2)
  json_list = []

  epi_fname = of_name + "epialign_res_" + runid
  seq_fname = of_name + "seqalign_res_" + runid
  epiScore_fname = of_name + "epi_scores_" + runid
  seqScore_fname = of_name + "seq_scores_" + runid
  out_name = of_name + "AlignResults_" + runid + ".txt"
  seq_stat = os.path.isfile(seq_fname)
  seq_bg = SeqBg(s, mu, alignMode)

  if seq_stat:
    fseq = open(seq_fname, "r")
  if alignMode == "enhancer":
    fepiScore = open(epiScore_fname, "r")
    if seq_stat:
      fseqScore = open(seqScore_fname, "r")

  with open(epi_fname, "r") as fepi, open(out_name, "w") as fout:
    i = 1
    line_seq = None
    line_seqScore = None
    while True:
      # Alignment results
      line_epi = fepi.readline().strip().split()
      if seq_stat:
        line_seq = fseq.readline().strip().split()
      # ALignment scores.
      if len(line_epi) == 0:
        break
      pair_name_raw = line_epi[0].split("_", 2)[-1]
      pair_name, one_num = pair_name_raw.split("$$$")

      try:
        #Initialize json_obj
        json_obj = InitJsonObj(i, pair_name, bed_dict1, bed_dict2, line_epi, line_seq, one_num)
        # Add region names.
        RegionName(json_obj, pair_name, intype1, intype2, alignMode)
        # The following steps are only for enhancer mode.
        if alignMode == "enhancer":
          query_len = json_obj["queryLength"]
          line_epiScore = [float(f) for f in fepiScore.readline().strip().split(",")[1:]]
          target_len = len(line_epiScore)
          line_epiScore = line_epiScore[0:target_len]
          if seq_stat:
            line_seqScore = [float(f) for f in fseqScore.readline().strip().split(",")[1:]]
            line_seqScore = line_seqScore[0:target_len]
          # Extract the two additional scores. Evaluate sequence similarity.

          SequenceEvaluation(json_obj, line_epi, line_seq, line_epiScore, line_seqScore, s, mu, seq_bg)

        # Write results to file.
        WriteFinalResult(json_obj, fout, alignMode)
        json_list.append(json_obj)
        i += 1
      except Exception as e:
        print >> sys.stderr, '[EpiAlignment]Error parsing item: ' + e.message

  json_dump_dict = {"data": json_list, "seqBg": seq_bg, "alignMode": alignMode, "searchRegionMode": searchRegionMode, "runid": runid}
  with open(of_name + runid + ".json", "w") as fjson:
    json.dump(json_dump_dict, fjson)
      


def Main():
  # Parse the json string passed by node js. 
  web_json = ParseJson()
  # Output folder name
  runid = web_json["runid"]
  allres_path = web_json["path"]
  out_folder = allres_path + "/tmp_" + web_json["runid"] + "/"

  sys.stdout.flush()
  # Move all uploaded files to the output folder. (this is done in server.js now)
  # MoveUploadFiles(out_folder, web_json["files"])

  # Generate input data
  # Parse peak files.
  ParsePeaks(out_folder, web_json, runid)
  # Create a pair of bed files.
  bed1, bed2, intype1, intype2 = CreateInputBeds(out_folder, web_json, runid)
  # Generate the input fastq-like file.
  BedToFa(bed1, bed2, out_folder, web_json["body"]["genomeAssembly"], runid)
  # Generate parameter file
  InputParas(out_folder, web_json["body"], runid)

  # Run EpiAlignment
  ExeEpiAlignment(web_json["body"]["alignMode"], web_json["body"]["searchRegionMode"], bed1, bed2, web_json["body"]["genomeAssembly"], out_folder, runid)
  # Parse the alignment results.
  ParseAlignResults(bed1, bed2, intype1, intype2, web_json["body"]["alignMode"], web_json["body"]["searchRegionMode"], \
    out_folder, runid, web_json["body"]["paras"], web_json["body"]["paramu"])








Main()
