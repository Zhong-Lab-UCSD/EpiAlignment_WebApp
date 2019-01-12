import requests
import json
import sys
import argparse
from time import time, sleep
import datetime

if sys.version_info[0] > 2:
  raise Exception("Please use Python 2 or convert the code to Python 3.")

EpiAligment_URL = "https://beta.epialign.ucsd.edu"

RUNNING_CODE = -1
# run info will be printed into the output file follow this order:
Run_info = [ ["runid", "alignMode", "subMode"],\
["queryGenomeAssembly", "queryInput", "queryPeak"],\
["targetGenomeAssembly", "clusters", "targetInput", "targetPeak"],\
["promoterUp", "promoterDown", "enhancerUp", "enhancerDown"],\
["seqWeight", "epiWeight", "paraS", "paraMu", "paraK", "piA", "piC", "piG", "piT", "pi1"] ]
# Results will be printed into the output file follow this order:
Out_list = ["index", "region_name1", "ensID1", "transID1", "region1",\
"region_name2", "ensID2", "transID2", "region2", "scoreE", "scoreS", "targetE", "targetS"]

def ParseArg():
  p = argparse.ArgumentParser(description = "The python client template for EpiAlignment.")
  p.add_argument("Input", type=str, nargs="?", help="Input file name. The input file is a sample sheet with details of jobs to be submitted to EpiAlignment.")
  p.add_argument("--public_data", action="store_true", help="If specified, get available ENCODE/public dataset IDs and descriptions in EpiAlignment.")
  p.add_argument("--find_gene_cluster", type=str, help="Search for gene clusters containing a gene name/Ensembl id. Please provide an ensembl id or a gene name/partial name.")
  if len(sys.argv) == 1:
    print >> sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args()


class SFObject:
  def __init__(self, fname):
    self.fin = open(fname, "r")
    self.header = self.fin.readline()

  def __iter__(self):
    return self

  def next(self):
    line = self.fin.next()
    return self.ParseFormLine(line)


  def ParseEssential(self, line, d_dict, f_list):
    # The speciesText represents inputs in textareas on the web page.
    d_dict["speciesText"] = ["", ""]
    
    # Parse shared elements.
    # Parse all parameters and add them to the data dictionary.   
    para_list = line[15:24]
    para_name_list = ["epiweight", "paras", "paramu", "parak", "piA", "piC", "piG", "piT", "pi1"]
    if "" in para_list or len(para_name_list) != len(para_list):
      print >> sys.stderr, "Please check the parameters."
      sys.exit(1)
    p_dict = zip(para_name_list, para_list)
    d_dict.update(p_dict)
    
    # Parse all peak files. If encode/public ids are provided, they will be added to the data dictionary. 
    # Otherwise, they will be added to the f_list and uploaded.
    # Note that the encode/public datasets will override user's peak files.
    peak_list = line[4:8]
    if peak_list[0] != "" and peak_list[1] != "":
      # encode/public data
      d_dict["encodeData"] = peak_list[0:2]
    elif peak_list[2] != "" and peak_list[3] != "":
      fpeak1 = peak_list[2]
      fpeak2 = peak_list[3]
      f_list.append(("speciesPeak1[]", open(fpeak1, "rb")))
      f_list.append(("speciesPeak2[]", open(fpeak2, "rb")))
    else:
      print >> sys.stderr, "Please provide paired peak files."
      sys.exit(1)
    
    # Add genome assemblies into the data dictionary.
    sp1 = line[2]
    sp2 = line[3]
    if sp1 != "" and sp2 != "":
      d_dict["genomeAssembly"] = [sp1, sp2]
    else:
      print >> sys.stderr, "Please provide genome assemblies."
      sys.exit(1)


  def ParseFormLine(self, line):
    '''
    Parse every line and create objects for the post request.
    Input: a list containing submission information, which is a line in the sample sheet.
    Return: a tuple contains a dict and a list. The first one contains data and the second one contains files.
    '''
    data_dict = {}
    file_list = []
    line = line.strip().split("\t")
    # Add elements including genome assemblies, peak files and parameters to the data dict.
    self.ParseEssential(line, data_dict, file_list)

    data_dict["alignMode"] = line[0]
    data_dict["searchRegionMode"] = line[1]
    finput1 = line[8]
    finput2 = line[9]
    cluster_name = line[10]
    promoterUp = line[11]
    promoterDown = line[12]
    enhancerUp = line[13]
    enhancerDown = line[14]
    # The first input file is required.
    if finput1 == "":
      print >> sys.stderr, "The query input file is required."
      sys.exit(1)
    file_list.append(("speciesInput1", open(finput1, "rb")))

    if line[1] == "genomeregion":
      # define target regions with an input file.
      if finput2 == "":
        print >> sys.stderr, "Please provide paired input files or change the searchRegionMode."
        sys.exit(1)
      file_list.append(("speciesInput2", open(finput2, "rb")))
      if line[0] == "promoter":
        data_dict["promoterUp"] = promoterUp
        data_dict["promoterDown"] = promoterDown

    else:
      if line[0] == "promoter" and line[1] == "genecluster":
        # search a gene cluster.
        if cluster_name == "":
          print >> sys.stderr, "Please provide a cluster name."
          sys.exit(1)
        data_dict["clusters"] = cluster_name
        data_dict["promoterUp"] = promoterUp
        data_dict["promoterDown"] = promoterDown

      elif line[0] == "enhancer" and line[1] == "homoregion":
        # define target regions using homologous regions.
        if enhancerUp == "" or enhancerDown == "":
          print >> sys.stderr, "Please provide enhancerUp and enhancerDown."
          sys.exit(1)
        data_dict["enhancerUp"] = enhancerUp
        data_dict["enhancerDown"] = enhancerDown

      else:
        print >> sys.stderr, "Please check your alignMode and searchRegionMode."
        sys.exit(1)
    return data_dict, file_list



class RequestsEpiAlign:
  def __init__(self, domain):
    self.domain = domain

  def get(self, path):
    return requests.get(self.domain + path)

  def post(self, path, **kwargs):
    # {runid: runid}.
    return requests.post(self.domain + path, **kwargs)

  def Post_sample(self, data, files):
    try:
      if len(files) != 0:
        r = self.post("/backend/form_upload", data = data, files = files)
      else:
        r = self.post("/backend/form_upload", data = data)
      r.raise_for_status()
    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
      print >> sys.stderr, "Fail to connect."
      return
    except requests.exceptions.HTTPError:
      print >> sys.stderr, "HTTP error."
      return
    else:
      runid_dict = json.loads(r.text)
      return runid_dict["runid"]

  def Get_sample(self, runid, gap=10, total_wait = 172800):
    '''
    Send get requests to the website repeatly until the task is done.

    Inputs:
    runid: the runid of current task, returned by the post request.
    gap: the time python will sleep between two get requests.
    total_wait: total waiting time.

    Return:
    A python dictionary. The results are stored in "data".
    '''
    time_start = time()
    time_waited = 0
    while time_waited < total_wait:
      try:
        r = self.get("/backend/results/" + runid)
        r.raise_for_status()
        json_dict = json.loads(r.text)
        if int(json_dict["status"]) == -1:
          # The task is still running. Wait for 'gap' time.
          print >> sys.stderr, "[" + str(datetime.datetime.now()) + "] Job " + runid + " is still running...\r",
          sleep(gap)
          time_waited += gap
        elif int(json_dict["status"]) == 0:
          # Finished successfully.
          print >> sys.stderr, ""
          print >> sys.stderr, "[" + str(datetime.datetime.now()) + "] Job " + runid + " has been finished successfully!"
          return json_dict
        else:
          # Finished with error.
          print >> sys.stderr, "Job finished with error. Error code: " + str(json_dict["status"])
          return

      except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
        print >> sys.stderr, "Fail to connect."
        return
      except requests.exceptions.HTTPError:
        print >> sys.stderr, "HTTP error."
        return

  def Get_cluster(self, partial_name):
    '''
    Search a gene name/ensemble id or a partial name to find
    gene clusters containing these names.

    Return: a dictionary.
    '''
    try:
      r = self.get("/backend/get_cluster/" + partial_name)
      r.raise_for_status()
      json_dict = json.loads(r.text)
      return json_dict

    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
      print >> sys.stderr, "Fail to connect."
      return
    except requests.exceptions.HTTPError:
      print >> sys.stderr, "HTTP error."
      return

  def Get_publicData(self):
    '''
    Get all ENCODE/public dataset IDs and descriptions in EpiAlignment.
    Return: two python dictionaries with encode datasets and public datasets information.
    '''
    try:
      r_encode = self.get("/assets/encodeData.json")
      r_public = self.get("/assets/publicData.json")
      r_encode.raise_for_status()
      r_public.raise_for_status()
      enc_dict = json.loads(r_encode.text)
      pub_dict = json.loads(r_public.text)
      return enc_dict, pub_dict
    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
      print >> sys.stderr, "Fail to connect."
      return
    except requests.exceptions.HTTPError:
      print >> sys.stderr, "HTTP error."
      return


def ParseOutput(json_dict, fout_name):
  '''
  Parse the output dictionary and print results into a tab-delimited file.
  Input: json_dict, the output dictionary of Get_sample(); fout_name: the output file name.
  '''

  with open(fout_name, "w") as fout:
    # Print task information into the file.
    # These lines all start with "#".
    for info_line in Run_info:
      line = "# "
      line += "  ".join([item + ": " + str(json_dict[item]) for item in info_line if item in json_dict])
      print >> fout, line

    # Print a header into the output file.
    print >> fout, "\t".join(["Index", "Query_region_name", "Query_gene", "Query_transcript", "Query_coordinate",\
     "Target_region_name", "Target_gene", "Target_transcript", "Target_coordinate", "EpiAlign_score",\
      "SeqOnly_score", "EpiAlign_target", "SeqOnly_target"])
    # Iterate region pairs and print results into the output file, one result per line.
    for res in json_dict["data"]:
      print >> fout, "\t".join([str(res[key]) for key in Out_list])


def PublicDataInfo(data_dict, fout):
  '''
  Extract dataset information from the output dictionaries of Get_publicData().
  Print the results into the output file, one dataset per line.

  Biosamples (i.e. tissue types and cell lines) are stored under the key
  matchedTissues in each dictionary as a list. Under each biosample, there are
  two species: human and mouse. Details of datasets of the two species are stored as
  dictionaries under each species name.
  '''
  # epiName_converter: a dictionary with readable names of epi marks and 
  epiName_converter = data_dict["filters"]
  for biosample in data_dict["matchedTissues"]:
    # biosample: a dictionary
    biosample_name = biosample[".label"]
    for species in biosample:
      if species[0] == ".":
        continue
      for dataset in biosample[species]["experiments"]:
        epiMark = epiName_converter[dataset["filterId"]]["target.label"]
        assay = epiName_converter[dataset["filterId"]]["assay_term_name"]
        if "biosampleId" in dataset:
          biosampleID = dataset["biosampleId"]
        else:
          biosampleID = "."
        peakId = dataset["id"]
        print >> fout, "\t".join([biosample_name, epiMark, assay, species, biosampleID, peakId])


def ParsePublicData(enc_dict, pub_dict):
  '''
  Parse details of the preset public datasets in EpiAlignment
  and print the information into a tab-delimited file named EpiAlign_publicData.txt

  Input: the output dictionaries of Get_publicData().
  '''
  fdata_name = "EpiAlign_publicData.txt"
  # column names
  dcol_names = ["biosample", "epiMark", "assay", "species", "biosampleId", "peakFileId"]
  with open(fdata_name, "w") as fout:
    print >> fout, "\t".join(dcol_names)
    PublicDataInfo(enc_dict, fout)
    PublicDataInfo(pub_dict, fout)


def ContainStr(gene, search_str):
  '''
  Check if the search_str matches the gene name or aliases.
  Print each gene.
  '''
  # gene symbols and ensembl ids will always be printed.
  contained_holder = "  * "
  uncontained_holder = "    "
  out_str = gene["symbol"] + " " + gene["ensemblId"]
  if search_str == gene["ensemblId"]:
    # if search string is an ensembl ID, it should always be fully matched to an ID in the cluster.
    print contained_holder + out_str
    return
  search_str = search_str.lower()
  if search_str in gene["symbol"].lower():
    print contained_holder + out_str
    return
  for altname in gene["aliases"]:
    if search_str in altname.lower():
      print contained_holder + out_str + " " + ",".join(gene["aliases"][1:])
      return
  print uncontained_holder + out_str


def PrintClusterGenes(cluster_list, search_str):
  '''
  Parse the cluster list and print genes in each cluster.
  '''
  for cluster in cluster_list:
    # cluster is a dictionary.
    # print cluster name
    print cluster["id"]
    for species in cluster["genesBySpecies"]:
      print species
      for gene in cluster["genesBySpecies"][species]:
        ContainStr(gene, search_str)
    print ""


def ParseClusterInfo(json_dict, search_str):
  '''
  Print found clusters and genes in each cluster.
  The partially or fully matched gene names/ids will be marked with an asterisk (*).

  Gene symbols and ensembl ids will always be printed out. If the query string matches
  an alias, all aliases of the gene will also be provided.
  '''
  print "Query string: %s"%(search_str)
  if json_dict["maxExceeded"]:
    print "The query string is too short. Please use a more specific query string to search for gene clusters."
    sys.exit(0)
  elif len(json_dict["fullMatchList"]) == 0 and len(json_dict["partialMatchList"]) == 0:
    print "The gene does not exist or there is no paralogues matched in selected species."
    sys.exit(0)
  else:
    if len(json_dict["fullMatchList"]) != 0:
      print ""
      print "Cluster(s) with fully-matched names:"
      PrintClusterGenes(json_dict["fullMatchList"], search_str)
    if len(json_dict["partialMatchList"]) != 0:
      print ""
      print "Cluster(s) with partially-matched names:"
      PrintClusterGenes(json_dict["partialMatchList"], search_str)
      


def Main():
  args = ParseArg()
  # Create a RequestsEpiAlign object to send http post and get requests.
  session = RequestsEpiAlign(EpiAligment_URL)
  if args.Input: 
    # Send jobs to Epialignment.
    fin_name = args.Input
    # Create a SFObject object to parse the input sample sheet. 
    # The return value SampleForm is an iterator. You may iterate it to 
    # get a data dictionary and a file list for each task (each row in your sample sheet)
    # for the post request.
    SampleForm = SFObject(fin_name)
    #
    # Start sending sample information to EpiAlignment.
    #
    # Please note that this program will not send another post request until the previous
    # task is done to prevent overloading the server. Please avoid sending several post 
    # requests simultaneously. 
    sample_index = 1
    for data, files in SampleForm:
      # Send a post requests to transfer data and files to EpiAlignment.
      # This function will return the runid of 
      runid = session.Post_sample(data, files)
      print runid
      # Get the result when the job is done. 
      # The function will return a python dictionary.
      # Use result_dict["data"] to access the list of results. 
      #
      # In this list, each element is a dictionary containing alignment scores and other information
      # of region pairs.
      result_dict = session.Get_sample(runid)
      # Print results into the output files. 
      fout_name = "alignResult_" + str(sample_index) + ".txt"
      ParseOutput(result_dict, fout_name)
      sample_index += 1
  elif args.public_data:
    # Get information of ENCODE/public datasets in EpiAlignment.
    # Details will be printed to a tab-delimited file named "EpiAlign_publicData.txt".
    #
    # You may use software such as excel to view the output file and select the datasets
    # you would like to use.
    encode_dict, public_dict = session.Get_publicData()
    ParsePublicData(encode_dict, public_dict)
  elif args.find_gene_cluster:
    # Search for gene clusters.
    cluster_dict = session.Get_cluster(args.find_gene_cluster)
    # Print found clusters and genes in each cluster.
    # 
    ParseClusterInfo(cluster_dict, args.find_gene_cluster)



if __name__ == '__main__':
	Main()