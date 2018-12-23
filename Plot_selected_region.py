import sys
import re
import os
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.vectors import FloatVector
import rpy2.robjects.packages as packages
from rpy2.robjects.vectors import IntVector
from rpy2.robjects import DataFrame
from rpy2.robjects.packages import importr
import numpy as np
rpy2.robjects.numpy2ri.activate()
grdevices = importr('grDevices')
rprint = ro.globalenv.get("print")

datasets = importr('datasets')
mtcars = packages.data(datasets).fetch('mtcars')['mtcars']

def Extract_name(fname, ind):
  with open(fname, "r") as fin:
    for line in fin:
      if line[0] == "#":
        continue
      line = line.strip().split("\t")
      if line[0] == "Index":
        continue
      if int(line[0]) == ind:
        xtitle = line[3].split('(')[0]
        coord = re.findall(r"[\w']+", xtitle)
        strand = line[3].split("(")[1].strip(")")
        return xtitle, int(coord[1]), int(coord[2]), strand


def Extract_selected_region(fname, ind):
  i = 1
  with open(fname, "r") as fin:
    for line in fin:
      if i == ind:
        line = line.strip().split(",")[1:]
        flag = 1
        return [float(f) for f in line]
      i += 1

  print >> sys.stderr, "No such image index."
  sys.exit(310)

def Plot_ScoreDist(s_list, ind, mode, of_name, runid):
  xtitle, start, stop, strand = Extract_name(of_name + "AlignResults_" + runid + ".txt", ind)
  tlen = stop - start
  if strand == "+":
    s_array = s_list[0:tlen]
  else:
    s_array = s_list[0:tlen][::-1]
  df = {"value": FloatVector(s_array), "coordinate": IntVector(list(range(start, stop)))}
  dataf = DataFrame(df)

  plot_name = of_name + "Image_" + str(ind) + "_" + mode + "_" + runid + ".png" 

  gp = ggplot2.ggplot(dataf)
  p = (gp + ggplot2.aes_string(x="coordinate", y='value') + ggplot2.geom_point(size = 0.1) +\
   ggplot2.labs(x=xtitle, y="Alignment scores") +\
   ggplot2.theme(text=ggplot2.element_text(size=14)))

  ro.r.ggsave(filename=plot_name, plot=p, width=150, height=60, unit='mm')
  #rprint(p)
  #grdevices.dev_off()


def Main():
  runid = sys.argv[1]
  ind = int(sys.argv[2])
  # mode can be "epi" or "seq"
  mode = sys.argv[3]
  #lines = sys.stdin.readlines()
  #runid = lines[0]

  lines = sys.stdin.readlines()
  json_dict = json.loads(lines[0])
  runid = json_dict["runid"]
  ind = int(json_dict["index"])
  mode = json_dict["mode"]

  out_folder = "tmp_" + runid + "/"
  if mode == "epi":
    fname = out_folder + "epi_scores_" + runid
    if not os.path.isfile(fname):
      sys.exit(0)
  elif mode == "seq":
    fname = out_folder + "seq_scores_" + runid
    if not os.path.isfile(fname):
      sys.exit(0)

  selected_list = Extract_selected_region(fname, ind)
  Plot_ScoreDist(selected_list, ind, mode, out_folder, runid)

Main()



