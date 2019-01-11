import subprocess
from subprocess import Popen, PIPE
from GeneAnno import *
import json
import shutil
import sys
import os
import re

ensembl_regexp = 'ENS[A-Z]+[0-9]{11}'

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


def CheckFileType(file1):
  with open(file1, "r") as fin:
    line = fin.readline().strip().split()
    if len(line) == 6:
      return "bed"
    elif len(line) == 1:
      return "name"
    else:
      print >> sys.stderr, "[EpiAlignment]Input files have to be bed6 files or genelists."
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

def FileOrTextarea(textarea_input, json_files, key, of_name, runid):
  '''
  Determine if input was pasted into the textarea or uploaded as a file.
  textarea_input: a string. Text in textarea.
  json_files: json "files" value.
  key: key value of the file.
  of_name: output folder name.
  return: a string and file type. If a file was uploaded, simply return file name. If data was pasted into the textarea,
  write the data into a new file.
  '''
  if key in json_files:
    # uploaded file
    fname = of_name + json_files[key][0]["filename"]
  elif textarea_input != "":
    # paste data
    fname = of_name + key + "_" + runid
    with open(fname, "w") as fsp:
      print >> fsp, textarea_input.rstrip("\n")
  else:
    # no data provided.
    return "", ""

  intype = CheckFileType(fname)
  return fname, intype

########################
## Mode 1: genelist   ##
########################
def Cons_transDict(gene_name, sp_name, of_name):
  if re.match(ensembl_regexp, gene_name):
    transDict = Construct_ensDict("Annotation/AnnotationFiles/" + sp_name + "_transcript.clean")
  else:
    transDict = Construct_nameDict("Annotation/AnnotationFiles/" + sp_name + "_transcript.clean")
  return transDict


def Cons_transList(input1, intype1, promoterUp, promoterDown, sp, of_name):
  trans_list1 = []
  with open(input1, "r") as fin1:
    if intype1 == "bed":
      trans_list1 = [line.strip().split() for line in fin1]
    elif intype1 == "name":
      i = 0
      for line in fin1:
        line = line.strip()
        if i == 0:
          transDict1 = Cons_transDict(line, sp, of_name)
          i += 1
        if line in transDict1:
          trans_list1 += PromoterMerge(line, transDict1, promoterUp, promoterDown)
        else:
          print >> sys.stderr, "[EpiAlignment]The gene " + line + " was not found in " + sp
  return trans_list1


def PairCutPromoter(input1, input2, intype1, intype2, promoterUp, promoterDown, genAssem, of_name):
  trans_list1 = Cons_transList(input1, intype1, promoterUp, promoterDown, genAssem[0], of_name)
  trans_list2 = Cons_transList(input2, intype2, promoterUp, promoterDown, genAssem[1], of_name)

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
  transDict2 = Cons_transDict(cluster_genes2[0], genAssem[1], of_name)
  trans_list2 = []
  for gene in cluster_genes2:
    if gene in transDict2:
      trans_list2 += PromoterMerge(gene, transDict2, promoterUp, promoterDown)
  
  fname2 = of_name + cluster_id + runid + "_2.bed"
  if input1 != "":
    # uploaded file
    trans_list1 = Cons_transList(input1, intype1, promoterUp, promoterDown, genAssem[0], of_name)
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
    transDict1 = Cons_transDict(cluster_genes1[0], genAssem[0], of_name)
    trans_list1 = []
    for gene in cluster_genes1:
      if gene in transDict1:
        trans_list1 += PromoterMerge(gene, transDic1, promoterUp, promoterDown)

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
      elif line[5] == "-":
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
    print >> sys.stderr, "[EpiAlignment]Failed to generate the input file. Exit code: " + str(exit_code)
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
  input1, intype1 = FileOrTextarea(json_dict["body"]["speciesText"][0], json_dict["files"], "speciesInput1", of_name, runid)
  if input1 == "" and (not json_dict["body"]["searchRegionMode"] == "genecluster"):
    print >> sys.stderr, "[EpiAlignment]No input regions provided."
    sys.exit(200)

  if searchMode == "genomeregion":
    # Mode 1: define search regions with bed files or gene lists.
    # Is input2 a file or a pasted text?
    input2, intype2 = FileOrTextarea(json_dict["body"]["speciesText"][1], json_dict["files"], "speciesInput2", of_name, runid)
    if alignMode == "enhancer":
      if CheckFileLength(input1, input2):
        return input1, input2, intype1, intype2
    else:
      bed1, bed2 = PairCutPromoter(input1, input2, intype1, intype2, promoterUp, promoterDown, genAssem, of_name)
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
      if StripDigits(genAssem[0]) == StripDigits(genAssem[1]):
        print >> sys.stderr, "[EpiAlignment]The two species must be different to use this mode."
        sys.exit(208)
      # Extend the input bed file1. extbed: extended bed file name.
      extbed = ExtendBed(input1, enhancerUp, enhancerDown)
      # LiftOver
      liftbed = LiftOver(extbed, genAssem)
      # Remove non-remappable regions. Return a pair of bed file names.
      cleanbed = RemoveNonlift(input1, liftbed)
      os.remove(extbed)
      os.remove(extbed + ".unlift")
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
  para_list = seq_pi_list + pi_list1 + weight_list
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


def ExeEpiAlignment(alignMode, searchRegionMode, of_name, runid):
  '''
  Execute EpiAlignment
  '''
  seq_stat = os.path.isfile(of_name + "parameters_seq_" + runid)
  cmd_list = ["python", "EpiAlignment.py", of_name + "Input_" + runid] +\
    ["-e", of_name + "parameters_" + runid] +\
    ["-p", "140"] +\
    ["-o", of_name + "epialign_res_" + runid]

  cmd_list_seq = ["python", "EpiAlignment.py", of_name + "Input_" + runid] +\
    ["-e", of_name + "parameters_seq_" + runid] +\
    ["-p", "140"] +\
    ["-o", of_name + "seqalign_res_" + runid]

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

def InitJsonObj(ind, pair_name, bed_dict1, bed_dict2, line_epi, line_seq = ""):
  '''
  Initialize a json object
  ifseq: If sequence-only alignment has been done. 
  '''
  json_obj = {"index":ind, "region1": ConcateBed(bed_dict1[pair_name]), "region2": ConcateBed(bed_dict2[pair_name]),\
  "scoreE": line_epi[1], "targetE": TargetRegion(bed_dict2[pair_name], line_epi),\
  "region_name1": ".", "region_name2": ".", "ensID1": ".", "ensID2": ".", "transID1":".", "transID2": "."}

  if line_seq == "":
    json_obj["scoreS"] = "."
    json_obj["targetS"] = "."
  else:
    json_obj["scoreS"] = line_seq[1]
    json_obj["targetS"] = TargetRegion(bed_dict2[pair_name], line_seq)

  return json_obj


def BedDict(fname):
  bed_dict = {}
  with open(fname, "r") as fin:
    for line in fin:
      line = line.strip().split()
      bed_dict[line[3]] = line[0:3] + [line[5]]
  return bed_dict


def TargetRegion(bed_list, res_line):
  if bed_list[3] == "+":
    start = int(bed_list[1]) + int(res_line[5])
    stop = int(bed_list[1]) + int(res_line[6])
  else:
    start = int(bed_list[2]) - int(res_line[6])
    stop = int(bed_list[2]) - int(res_line[5])
  return bed_list[0] + ":" + str(start) + "-" + str(stop) + "(" + bed_list[3] + ")"


def ConcateBed(coor_list):
  return coor_list[0] + ":" + str(coor_list[1]) + "-" + str(coor_list[2]) + "(" + coor_list[3] + ")"


def WriteFinalResult(json_obj, fout, alignMode):
  if json_obj["index"] == 1:
    print >> fout, "\t".join(["Index", "Query_region_name", "Query_gene", "Query_transcript", "Query_coordinate",\
     "Target_region_name", "Target_gene", "Target_transcript", "Target_coordinate", "EpiAlign_score",\
      "SeqOnly_score", "EpiAlign_target", "SeqOnly_target"])

  print >> fout, "\t".join([str(f) for f in [json_obj["index"],\
    json_obj["region_name1"], json_obj["ensID1"], json_obj["transID1"], json_obj["region1"],\
    json_obj["region_name2"], json_obj["ensID2"], json_obj["transID2"], json_obj["region2"],\
    json_obj["scoreE"], json_obj["scoreS"], json_obj["targetE"], json_obj["targetS"] ] ])

def ParseAlignResults(bed1, bed2, intype1, intype2, alignMode, searchRegionMode, of_name, runid):
  '''
  Parse alignment results.
  bed1, bed2: the two bed files used for generating input file.
  intype: "bed", "name" or an empty string.
  of_name: output folder.
  runid: runid.
  return: None. This function will write a json object to a file.
  '''
  bed_dict1 = BedDict(bed1)
  bed_dict2 = BedDict(bed2)
  json_list = []

  epi_fname = of_name + "epialign_res_" + runid
  seq_fname = of_name + "seqalign_res_" + runid
  out_name = of_name + "AlignResults_" + runid + ".txt"
  seq_stat = os.path.isfile(seq_fname)

  if seq_stat:
    fseq = open(seq_fname, "r")

  with open(epi_fname, "r") as fepi, open(out_name, "w") as fout:
    i = 1
    while True:
      line_epi = fepi.readline().strip().split()
      if seq_stat:
        line_seq = fseq.readline().strip().split()
      if len(line_epi) == 0:
        break
      pair_name = line_epi[0].split("_", 2)[-1]

      #Initialize json_obj
      if seq_stat:
        json_obj = InitJsonObj(i, pair_name, bed_dict1, bed_dict2, line_epi, line_seq)
      else:
        json_obj = InitJsonObj(i, pair_name, bed_dict1, bed_dict2, line_epi)

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

      WriteFinalResult(json_obj, fout, alignMode)
      json_list.append(json_obj)
      i += 1

  json_dump_dict = {"data": json_list, "alignMode": alignMode, "searchRegionMode": searchRegionMode, "runid": runid}
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
  # Move all uploaded files to the output folder.
  MoveUploadFiles(out_folder, web_json["files"])

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
  ExeEpiAlignment(web_json["body"]["alignMode"], web_json["body"]["searchRegionMode"], out_folder, runid)
  # Parse the alignment results.
  ParseAlignResults(bed1, bed2, intype1, intype2, web_json["body"]["alignMode"], web_json["body"]["searchRegionMode"], out_folder, runid)








Main()