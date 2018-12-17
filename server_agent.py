import subprocess
import json
import shutil
import sys

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


def ParsePeaks(of_name, json_files, runid):
  '''
  Creak peak files.
  '''
  with open(of_name + "peaks_" + runid) as fpeak:
    print >> fpeak, "@species1"
    print >> fpeak, of_name + json_files["speciesPeak[]"][0]["filename"]
    print >> fpeak, "@species2"
    print >> fpeak, of_name + json_files["speciesPeak[]"][1]["filename"]


def CreateInputBeds(of_name, json_dict):
  '''
  Create a pair of bed files.
  of_name: output folder name.
  json_dict: the json dictionary.
  return: names of the two bed files.
  '''

  input1 = json_dict["files"]["speciesInput[]"][0]["filename"]

  if web_json["body"]["searchRegionMode"] == "genomregion":
    # Mode 1: define search regions with bed files or gene lists.
    input2 = json_dict["files"]["speciesInput[]"][1]["filename"]
    if CheckFileLength(ofname + input1, of_name + input2):
      intype1 = CheckFileType(input1)
      intype2 = CheckFileType(input2)
      if intype1 == "bed" and intype2 == "bed":
        # Two input files are all bed6:
        bed1 = input1
        bed2 = input2
      else:
        # at least one file is a gene list. Cut the promoter regions.
        pass




  else:    
    if web_json["body"]["searchRegionMode"] == "genetype" and web_json["body"]["alignMode"] == "promoter":
      # Mode 2: search the promoter regions of a specific type of gene.
      pass
    elif web_json["body"]["searchRegionMode"] == "genecluster" and web_json["body"]["alignMode"] == "promoter":
      # Mode 3: search a specific gene cluster.
      pass

    elif web_json["body"]["searchRegionMode"] == "homoregion" and web_json["body"]["alignMode"] == "enhancer":
      # Mode 4 (enhancer mode 2): use homologous regions. 
      # species must be different!
      pass

  return bed1, bed2


def BedToFa(bed1, bed2, out_folder, sp_list, epiName_list, runid):
  '''
  Run InputToFastq_bed2.py to convert the bed file pair to a fastq-like file.
  '''
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
  with open(of_name + "parameters_" + runid) as fpara:
    print >> fpara, json_body["paras"]
    print >> fpara, json_body["paramu"]
    print >> fpara, "\n".join(json_body["parak"].split(","))
    print >> fpara, "A:" + json_body["piA"] + "\t" + "C:" + json_body["piC"] + "\t" +\
    "G:" + json_body["piG"] + "\t" + "T:" + json_body["piT"]
    kappa_list1 = [int(k) for k in json_body["pi1"].split(",")]
    kappa_list0 = [1 - k for k in kappa_list1]
    for p0, p1 in zip(kappa_list0, kappa_list1):
      print >> fpara, "0:" + str(p0) + "\t" + "1:" + str(p1)
    weights = "\t".join(json_body["epiweight"].split(","))
    print >> fpara, json_body["seqweight"] + "\t" + weights

def ExeEpiAlignment(of_name, runid):
  '''
  Execute EpiAlignment
  '''
  cmd_list = ["python", "EpiAlignment.py", out_folder + "Input_" + runid] +\
  ["-e", out_folder + "parameters_" + runid] +\
  ["-p", "140"] +\
  ["-o", out_folder + "epialign_result_" + runid]

  exit_code = subprocess.call(cmd_list)

  if exit_code != 0:
    print >> sys.stderr, "Failed to align regions. Exit code: " + str(exit_code)
    sys.exit(exit_code)


def Main():
  # Parse the json string passed by node js. 
  web_json = ParseJson()
  # Output folder name
  out_folder = web_json["runid"] + "/"
  # Move all uploaded files to the output folder.
  print web_json["files"]
  MoveUploadFiles(out_folder, web_json["files"])

  # Generate input data
  # Parse peak files.
  ParsePeaks(out_folder, web_json["files"], runid)
  # Create a pair of bed files.
  bed1, bed2 = CreateInputBeds(out_folder, web_json)
  # Generate the input fastq-like file.
  BedToFa(bed1, bed2, out_folder, web_json["body"]["genomeAssembly"], web_json["body"]["epiName"], web_json["body"]["runid"])
  # Generate parameter file
  InputParas(out_folder, web_json["body"], runid)

  # Run EpiAlignment
  ExeEpiAlignment(out_folder, runid)






Main()