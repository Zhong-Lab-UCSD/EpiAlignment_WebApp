import sys
import json
import re
import os
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects.packages as packages
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import DataFrame
from rpy2.robjects.packages import importr
import numpy as np
from math import *
rpy2.robjects.numpy2ri.activate()
gridExtra = importr("gridExtra")

GRAY_COLOR = "gray60"
BLACK_COLOR = "black"

def Extract_name(fname, ind):
  with open(fname, "r") as fin:
    for line in fin:
      if line[0] == "#":
        continue
      line = line.strip().split("\t")
      if line[0] == "Index":
        continue
      if int(line[0]) == ind:
        xtitle = line[8].split('(')[0]
        coord = re.findall(r"[\w']+", xtitle)
        strand = line[8].split("(")[1].strip(")")
        return xtitle, int(coord[1]), int(coord[2]), strand


def Extract_selected_region(fname, ind, tlen):
  '''
  tlen: length of the target region.
  '''
  i = 1
  with open(fname, "r") as fin:
    for line in fin:
      if i == ind:
        line = line.strip().split(",")[1:]
        query_len = len(line) - tlen
        norm_factor = 1000.0 / query_len
        return [float(f) * norm_factor for f in line]
      i += 1

  print >> sys.stderr, "No such image index."
  sys.exit(310)

def averageList(list):
  return float(sum(list)) / len(list)

def centerOfList(list):
  return list[int((len(list) - 1) / 2)]

def reduceList(list, reduceBy, func):
  result = []
  for i in xrange(int(ceil((len(list) - 1) / reduceBy)) + 1):
    result.append(func(list[(i * reduceBy):((i + 1) * reduceBy)]))
  return result

def Plot_ScoreDist(epi_list, seq_list, ind, of_name, runid, xtitle, start, stop, strand, tlen, imgFormat, reduceBy):
  seq_stat = 1 if seq_list != "" else 0
  reduceBy = 1 if int(reduceBy) < 1 else int(reduceBy)
  
  qlen = len(epi_list) - tlen
  # concatenate the last bit correctly
  epi_array = epi_list[0:tlen] + epi_list[tlen:][::-1]
  seq_array = seq_list[0:tlen] + seq_list[tlen:][::-1]

  # correct strand so the values match genomic coordinates
  if strand == "-":
    epi_array = epi_array[::-1]
    if seq_stat:
      seq_array = seq_array[::-1]
  
  coorRange = list(range(int(start - qlen / 2), int(ceil(stop + qlen / 2.0))))

  # Separate partial overlap (first and last `qlen` length) with full overlap
  
  if not seq_stat:
    df = {
      "value": FloatVector(reduceList(epi_array, reduceBy, averageList)),
      "coordinate": IntVector(reduceList(coorRange, reduceBy, averageList)),
      "mode": StrVector(
        reduceList(["EpiAlignment"] * (tlen + qlen), reduceBy, centerOfList)
      )
    }
  else:
    df = {
      "value": FloatVector(
        reduceList(epi_array, reduceBy, averageList) +
        reduceList(seq_array, reduceBy, averageList)),
      "coordinate": IntVector(reduceList(coorRange, reduceBy, averageList) * 2),
      "mode": StrVector(
        reduceList(["EpiAlignment"] * (tlen + qlen), reduceBy, centerOfList) +
        reduceList(["Sequence-only"] * (tlen + qlen), reduceBy, centerOfList))
    }
 
  if qlen < tlen:
    colorVec = [GRAY_COLOR] * qlen + [BLACK_COLOR] * (tlen - qlen) +\
      [GRAY_COLOR] * qlen
  else:
    colorVec = [GRAY_COLOR] * (tlen + qlen)
  
  colorVec = StrVector(
    reduceList(colorVec, reduceBy, centerOfList) *
    (2 if seq_stat else 1)
  )

  dataf = DataFrame(df)

  plot_name = of_name + "Image_" + str(ind) + "_" + runid + "." +\
    ("pdf" if imgFormat == "pdf" else "png")

  gp = ggplot2.ggplot(dataf)
  p = (gp + ggplot2.aes_string(x="coordinate", y='value') + \
    ggplot2.geom_point(size = 0.1, color = colorVec, shape = ".") +\
    ggplot2.labs(x=xtitle, y="Alignment scores") +\
    ggplot2.facet_grid(ro.Formula('mode ~ .'), scales = "free_y") +\
    ggplot2.theme(text=ggplot2.element_text(size=14)))

  if seq_stat:
    ro.r.ggsave(filename=plot_name,
                device="pdf" if imgFormat == "pdf" else "png",
                plot=p, width=150, height=120, unit='mm')
  else:
    ro.r.ggsave(filename=plot_name,
                device="pdf" if imgFormat == "pdf" else "png",
                plot=p, width=150, height=60, unit='mm')

def Main():

  lines = sys.stdin.readlines()
  json_dict = json.loads(lines[0])
  #runid = sys.argv[1]
  #ind = int(sys.argv[2])
  allpath_res = json_dict["path"]
  runid = json_dict["runid"]
  ind = int(json_dict["index"])
  imgFormat = json_dict["format"]
  reduceBy = json_dict["reduceBy"]
 

  out_folder = allpath_res + "/tmp_" + runid + "/"

  # Extract 
  xtitle, start, stop, strand = Extract_name(out_folder + "AlignResults_" + runid + ".txt", ind)
  tlen = stop - start
 

  fename = out_folder + "epi_scores_" + runid
  if not os.path.isfile(fename):
    sys.exit(0)
  selected_list_epi = Extract_selected_region(fename, ind, tlen)

  fsname = out_folder + "seq_scores_" + runid
  if not os.path.isfile(fsname):
    selected_list_seq = ""
  else:
    selected_list_seq = Extract_selected_region(fsname, ind, tlen)

  Plot_ScoreDist(selected_list_epi, selected_list_seq, ind, out_folder,
    runid, xtitle, start, stop, strand, tlen, imgFormat, reduceBy)

Main()



