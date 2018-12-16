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
      print ufile["path"]
      print dest_name + ufile["filename"]
      shutil.move(ufile["path"], dest_name + ufile["filename"])



def Main():
  # Parse the json string passed by node js. 
  web_json = ParseJson()
  # Output folder name
  out_folder = web_json["runid"]
  # Move all uploaded files to the output folder.
  print web_json["files"]
  MoveUploadFiles(out_folder, web_json["files"])



  if web_json["body"]["alignMode"] == "promoter":
    # Promoter mode
    
    pass
  elif web_json["body"]["alignMode"] == "enhancer":
    pass


Main()