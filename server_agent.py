import subprocess
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
  print "Linked to python..."
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
      print >> sys.stderr, "The two input files have different numbers of regions."
      sys.exit(202)


def CheckFileType(file1):
  with open(file1, "r") as fin:
    line = fin.readline().strip().split("\t")
    if len(line) == 6:
      return "bed"
    elif len(line) == 1:
      return "name"
    else:
      print >> sys.stderr, "Input files have to be bed6 files or genelists."
      sys.exit(201)

def StripDigits(qstr):
  '''
  Remove all digits in qstr.
  '''
  return ''.join([i for i in qstr if not i.isdigit()])


def ParsePeaks(of_name, json_files, runid):
  '''
  Creak peak files.
  '''
  if "speciesPeak1[]" not in json_files or "speciesPeak2[]" not in json_files :
    print >> sys.stderr, "No peak files found!"
    sys.exit(207)

  with open(of_name + "peaks_" + runid, "w") as fpeak:
    print >> fpeak, "@species1"
    print >> fpeak, of_name + json_files["speciesPeak1[]"][0]["filename"]
    print >> fpeak, "@species2"
    print >> fpeak, of_name + json_files["speciesPeak2[]"][0]["filename"]

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
    print >> sys.stderr, "No input regions provided."
    sys.exit(200)

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


def PairCutPromoter(input1, input2, intype1, intype2, promoterUp, promoterDown, genAssem, of_name):
  if intype1 == "bed" and intype2 == "bed":
    return input1, input2

  with open(input1, "r") as fin1, open(input2, "r") as fin2, open(input1 + ".bed", "w") as fout1, open(input2 + ".bed", "w") as fout2:
    i = 0
    while True:
      line1 = fin1.readline().strip()
      line2 = fin2.readline().strip()
      if line1 == "" or line2 == "":
        break
      if i == 0:
        # first line
        if intype1 == "name":
          transDict1 = Cons_transDict(line1, genAssem[0], of_name)
        if intype2 == "name":
          transDict2 = Cons_transDict(line2, genAssem[1], of_name)
        i += 1

      if intype1 == "name":
        if line1 not in transDict1:
          # No such gene
          continue
        trans_list1 = [PromoterBed(x, promoterUp, promoterDown) for x in transDict1[line1]]
      else:
        trans_list1 = [line1.split("\t")]

      if intype2 == "name":
        if line2 not in transDict2:
          continue
        trans_list2 = [PromoterBed(x, promoterUp, promoterDown) for x in transDict2[line2]]
      else:
        trans_list2 = [line2.split("\t")]

      for region1 in trans_list1:
        for region2 in trans_list2:
          region_name = region1[3] + "[===]" + region2[3]
          print >> fout1, "\t".join(region1[0:3] + [region_name] + region1[4:])
          print >> fout2, "\t".join(region2[0:3] + [region_name] + region2[4:])

  return input1 + ".bed", input2 + ".bed"


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
      line = line.strip().split("\t")
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
    print >> sys.stderr, "Failed to generate the input file. Exit code: " + str(exit_code)
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
          print >> sys.stderr, "None of the input regions were remappable."
          sys.exit(209)
        break
      i += 1
      line2 = line2.split("\t")

      while True:
        line1 = fin1.readline().strip().split("\t")
        if line1[3] == line2[3]:
          print >> fout, "\t".join(line1) 
          break
  return input_bed + ".clean"

##################
## Create beds  ##
##################

def CreateInputBeds(of_name, json_dict, runid):
  '''
  Create a pair of bed files.
  of_name: output folder name.
  json_dict: the json dictionary.
  return: names of the two bed files.
  '''
  # Is input1 a file or a pasted text?
  input1, intype1 = FileOrTextarea(json_dict["body"]["speciesText"][0], json_dict["files"], "speciesInput1", of_name, runid)

  if json_dict["body"]["searchRegionMode"] == "genomeregion":
    # Mode 1: define search regions with bed files or gene lists.
    # Is input2 a file or a pasted text?
    input2, intype2 = FileOrTextarea(json_dict["body"]["speciesText"][1], json_dict["files"], "speciesInput2", of_name, runid)
    if CheckFileLength(input1, input2):
      bed1, bed2 = PairCutPromoter(input1, input2, intype1, intype2, int(json_dict["body"]["promoterUp"]), int(json_dict["body"]["promoterDown"]), json_dict["body"]["genomeAssembly"], of_name)
    return bed1, bed2, intype1, intype2

  else:    
    if json_dict["body"]["searchRegionMode"] == "genetype" and json_dict["body"]["alignMode"] == "promoter":
      # Mode 2: search the promoter regions of a specific type of gene.
      pass

    elif json_dict["body"]["searchRegionMode"] == "genecluster" and json_dict["body"]["alignMode"] == "promoter":
      # Mode 3: search a specific gene cluster.
      pass

    elif json_dict["body"]["searchRegionMode"] == "homoregion" and json_dict["body"]["alignMode"] == "enhancer":
      # Mode 4 (enhancer mode 2): use homologous regions. 
      # species must be different!
      if StripDigits(json_dict["body"]["genomeAssembly"][0]) == StripDigits(json_dict["body"]["genomeAssembly"][1]):
        print >> sys.stderr, "The two species must be different to use this mode."
        sys.exit(208)
      # Extend the input bed file1. extbed: extended bed file name.
      extbed = ExtendBed(input1, int(json_dict["body"]["enhancerUp"]), int(json_dict["body"]["enhancerDown"]))
      # LiftOver
      liftbed = LiftOver(extbed, json_dict["body"]["genomeAssembly"])
      # Remove non-remappable regions. Return a pair of bed file names.
      cleanbed = RemoveNonlift(input1, liftbed)
      os.remove(extbed)
      os.remove(extbed + ".unlift")
      bed1 = cleanbed
      bed2 = liftbed

    return bed1, bed2, intype1, ""


def BedToFa(bed1, bed2, out_folder, sp_list, epiName, runid):
  '''
  Run InputToFastq_bed2.py to convert the bed file pair to a fastq-like file.
  '''
  epiName_list = epiName.split(",")
  cmd_list = ["python", "InputToFastq_bed2.py", bed1, bed2, "-s"] + sp_list +\
  ["--bg", out_folder + "peaks_" + runid] +\
  ["--histone"] + epiName_list +\
  ["-p", "20"] +\
  ["-o", out_folder + "Input_" + runid]

  exit_code = subprocess.call(cmd_list)

  if exit_code != 0:
    print >> sys.stderr, "Failed to generate the input file. Exit code: " + str(exit_code)
    sys.exit(exit_code)


def InputParas(of_name, json_body, runid):
  '''
  Create input parameter file for EpiAlignment.
  '''
  # Check if parameters are not positive.
  parak_list = [float(x) for x in json_body["parak"].split(",")]
  if float(json_body["paras"]) <= 0 or float(json_body["paramu"]) <= 0 or min(parak_list) <= 0:
    print >> sys.stderr, "Parameters should be positive values."
    sys.exit(206)

  seq_pi_list = [float(json_body["piA"]), float(json_body["piC"]), float(json_body["piG"]), float(json_body["piT"])]
  pi_list1 = [float(k) for k in json_body["pi1"].split(",")]
  weight_list = [float(w) for w in json_body["epiweight"].split(",")]
  para_list = seq_pi_list + pi_list1 + weight_list
  if min(para_list) <= 0 or max(para_list) >= 1:
    print >> sys.stderr, "Equilibrium probabilities (pi) must be values between 0 and 1."
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
    print >> fpara, json_body["seqweight"] + "\t" + weights

  # If epi weight is not 0, create another parameter file for sequence-only alignment.
  if float(json_body["seqweight"]) != 1:
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
  if alignMode == "promoter":
    cmd_list = ["python", "EpiAlignment.py", of_name + "Input_" + runid] +\
    ["-e", of_name + "parameters_" + runid] +\
    ["-p", "140"] +\
    ["-o", of_name + "epialign_res_" + runid]

    exit_code = subprocess.call(cmd_list)
    if exit_code != 0:
      print >> sys.stderr, "Failed to align regions. Exit code: " + str(exit_code)
      sys.exit(exit_code)

    if os.path.isfile(of_name + "parameters_seq_" + runid):
      cmd_list_seq = ["python", "EpiAlignment.py", of_name + "Input_" + runid] +\
      ["-e", of_name + "parameters_seq_" + runid] +\
      ["-p", "140"] +\
      ["-o", of_name + "seqalign_res_" + runid]

      exit_code = subprocess.call(cmd_list_seq)
      if exit_code != 0:
        print >> sys.stderr, "Failed to align regions. Exit code: " + str(exit_code)
        sys.exit(exit_code)

  elif alignMode == "enhancer":
    cmd_list = ["python", "EpiAlignment.py", of_name + "Input_" + runid] +\
    ["-e", of_name + "parameters_" + runid] +\
    ["-p", "140"] +\
    ["-o", of_name + "epialign_res_" + runid] +\
    ["-O", of_name + "epi_scores_" + runid]

    exit_code = subprocess.call(cmd_list)
    if exit_code != 0:
      print >> sys.stderr, "Failed to align regions. Exit code: " + str(exit_code)
      sys.exit(exit_code)

    if os.path.isfile(of_name + "parameters_seq_" + runid):
      cmd_list_seq = ["python", "EpiAlignment.py", of_name + "Input_" + runid] +\
      ["-e", of_name + "parameters_seq_" + runid] +\
      ["-p", "140"] +\
      ["-o", of_name + "seqalign_res_" + runid] +\
      ["-O", of_name + "seq_scores_" + runid]

      exit_code = subprocess.call(cmd_list_seq)
      if exit_code != 0:
        print >> sys.stderr, "Failed to align regions. Exit code: " + str(exit_code)
        sys.exit(exit_code)

###################
## Parse results ##
###################

def BedDict(fname):
  bed_dict = {}
  with open(fname, "r") as fin:
    for line in fin:
      line = line.strip().split("\t")
      bed_dict[line[3]] = line[0:3]
  return bed_dict


def TargetRegion(bed_list, res_line):
  start = int(bed_list[1]) + int(res_line[5])
  stop = int(bed_list[1]) + int(res_line[6])
  return bed_list[0] + ":" + str(start) + "-" + str(stop)


def ConcateBed(coor_list):
  return coor_list[0] + ":" + str(coor_list[1]) + "-" + str(coor_list[2])


def WriteFinalResult(json_obj, fout, alignMode):
  if alignMode == "enhancer":
    if json_obj["index"] == 1:
      print >> fout, "\t".join(["Index", "Region_name", "Query_coordinate", "Search_coordinate", "EpiAlign_score",\
        "SeqOnly_score", "EpiAlign_target", "SeqOnly_target"])

    print >> fout, "\t".join([str(f) for f in [json_obj["index"], json_obj["region_name"],\
      json_obj["region1"], json_obj["region2"],\
      json_obj["scoreE"], json_obj["scoreS"], json_obj["targetE"], json_obj["targetS"] ] ])
  else:
    if "ensID1" not in json_obj and "ensID2" not in json_obj:
      if json_obj["index"] == 1:
        print >> fout, "\t".join(["Index", "Region_name", "Query_coordinate", "Search_coordinate", "EpiAlign_score",\
          "SeqOnly_score"])

      print >> fout, "\t".join([str(f) for f in [json_obj["index"], json_obj["region_name"],\
        json_obj["region1"], json_obj["region2"],\
        json_obj["scoreE"], json_obj["scoreS"] ] ])

    elif "ensID1" not in json_obj and "ensID2" in json_obj:
      if json_obj["index"] == 1:
        print >> fout, "\t".join(["Index", "Query_region", "Query_coordinate", "Search_gene", "Search_transcript",\
          "Search_coordinate", "EpiAlign_score", "SeqOnly_score"])

      print >> fout, "\t".join([str(f) for f in [json_obj["index"], json_obj["region_name1"],\
        json_obj["region1"], json_obj["ensID2"], json_obj["transID2"], json_obj["region2"],\
        json_obj["scoreE"], json_obj["scoreS"] ] ])

    elif "ensID1" in json_obj and "ensID2" not in json_obj:
      if json_obj["index"] == 1:
        print >> fout, "\t".join(["Index", "Query_gene", "Query_transcript", "Query_coordinate", "Search_region",\
          "Search_coordinate", "EpiAlign_score", "SeqOnly_score"])

      print >> fout, "\t".join([str(f) for f in [json_obj["index"], json_obj["ensID1"], json_obj["transID1"],\
        json_obj["region1"], json_obj["region_name2"], json_obj["region2"],\
        json_obj["scoreE"], json_obj["scoreS"] ] ])

    else:
      if json_obj["index"] == 1:
        print >> fout, "\t".join(["Index", "Query_gene", "Query_transcript", "Query_coordinate", "Search_gene",\
          "Search_transcript", "Search_coordinate", "EpiAlign_score", "SeqOnly_score"])

      print >> fout, "\t".join([str(f) for f in [json_obj["index"], json_obj["ensID1"], json_obj["transID1"],\
        json_obj["region1"], json_obj["ensID2"], json_obj["transID2"], json_obj["region2"],\
        json_obj["scoreE"], json_obj["scoreS"] ] ])



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
  with open(of_name + "epialign_res_" + runid, "r") as fepi, open(of_name + "seqalign_res_" + runid, "r") as fseq,\
  open(of_name + "AlignResults_" + runid + ".txt", "w") as fout:
    # index, region1, region2, scoreS, scoreE.
    i = 1
    while True:
      line_epi = fepi.readline().strip().split("\t")
      line_seq = fseq.readline().strip().split("\t")
      if line_epi[0] == "":
        break
      pair_name = line_epi[0].split("_", 2)[-1]
      json_obj = {"index":i, "region1": ConcateBed(bed_dict1[pair_name]), "region2": ConcateBed(bed_dict2[pair_name]),\
       "scoreE": line_epi[1], "scoreS": line_seq[1],\
       "targetE": TargetRegion(bed_dict2[pair_name], line_epi), "targetS": TargetRegion(bed_dict2[pair_name], line_seq)}
      if (alignMode == "enhancer") or (intype1 == "bed" and intype2 == "bed"):
        # a pair of bed.
        json_obj["region_name"] = pair_name

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

  with open(of_name + runid + ".json", "w") as fjson:
    json.dump({"data": json_list}, fjson)
      


def Main():
  # Parse the json string passed by node js. 
  web_json = ParseJson()
  # Output folder name
  runid = web_json["runid"].split("_")[1]
  out_folder = web_json["runid"] + "/"
  print runid
  # Move all uploaded files to the output folder.
  MoveUploadFiles(out_folder, web_json["files"])

  # Generate input data
  # Parse peak files.
  ParsePeaks(out_folder, web_json["files"], runid)
  # Create a pair of bed files.
  bed1, bed2, intype1, intype2 = CreateInputBeds(out_folder, web_json, runid)
  # Generate the input fastq-like file.
  BedToFa(bed1, bed2, out_folder, web_json["body"]["genomeAssembly"], web_json["body"]["epiName"], runid)
  # Generate parameter file
  InputParas(out_folder, web_json["body"], runid)

  # Run EpiAlignment
  ExeEpiAlignment(web_json["body"]["alignMode"], web_json["body"]["searchRegionMode"], out_folder, runid)
  # Parse the alignment results.
  ParseAlignResults(bed1, bed2, intype1, intype2, web_json["body"]["alignMode"], web_json["body"]["searchRegionMode"], out_folder, runid)








Main()