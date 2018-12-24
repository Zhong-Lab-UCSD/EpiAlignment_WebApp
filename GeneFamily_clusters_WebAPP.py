import sys, argparse, os
from xplib.Annotation import Bed
from xplib import DBI

class OrthoRegion:
    def __init__(self, h_list):
        self.id1 = Unique([x[0] for x in h_list])
        self.id2 = Unique([x[2] for x in h_list])

class Paralog_Cluster:
    def __init__(self, id_list, name_list):
        self.id = id_list
        self.name = name_list
    def AddGene(self, newid, newname):
        self.id.append(newid)
        self.name.append(newname)


def ParseArg():
    p=argparse.ArgumentParser(description="Find homologous gene clusters.")
    p.add_argument("-s", "--species", nargs=2, default=["human", "mouse"], help="Species names.")
    p.add_argument("-r", "--paralog", nargs=2, type=str, help="Paralog files of the two species.")
    p.add_argument("-t", "--ortho", nargs=2, type=str, help="Homolog files of the two sepcies. Files should have header and four columns: GeneID1, GeneName1, GeneID2, GeneName2, Homology type.")
    p.add_argument("-o", "--output", type=str, help="output file name")
    if len(sys.argv) == 1:
        print >>sys.stderr, p.print_help()
        sys.exit(0)
    return p.parse_args() 

def Unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def Unique2(seq1, seq2):
    seen = set()
    seen_add = seen.add
    return [y for (x,y) in zip(seq1, seq2) if not (x in seen or seen_add(x))] 


def ParseOrtho(fname):
    tmp_dict = {}
    ortho_dict = {}
    with open(fname, "r") as fin:
        fin.readline()
        for line in fin:
            line = line.strip().split("\t")
            if line[2] != "" and line[5] != "NA":
                if line[0] not in tmp_dict:
                    tmp_dict[line[0]] = []
                tmp_dict[line[0]].append(line)

    for key in tmp_dict:
        ortho_dict[key] = OrthoRegion(tmp_dict[key])

    return ortho_dict

def ParseParalog(fname):
    '''
    para_dict: key: set name, value: a cluster
    gene_para_dict: key: gene ens_id, value: set name.
    '''
    with open(fname, "r") as fin:
        fin.readline()
        i = 1
        para_dict = {}
        gene_para_dict = {}
        for line in fin:
            line = line.strip().split("\t")
            if len(line) == 2:
                # No paralogous gene
                set_name = "Set" + str(i)
                para_dict[set_name] = Paralog_Cluster([line[0]], [line[1]])
                gene_para_dict[line[0]] = set_name
                i += 1
            elif line[0] not in gene_para_dict and line[2] not in gene_para_dict:
                # New set.
                set_name = "Set" + str(i)
                para_dict[set_name] = Paralog_Cluster([line[0], line[2]], [line[1], line[3]])
                gene_para_dict[line[0]] = set_name
                gene_para_dict[line[2]] = set_name
                i += 1
            elif line[0] in gene_para_dict and line[2] in gene_para_dict:
                # duplicated line with opposite order.
                continue
            else:
                if line[2] in gene_para_dict:
                    set_name = gene_para_dict[line[2]]
                    para_dict[set_name].AddGene(line[0], line[1])
                    gene_para_dict[line[0]] = set_name
                elif line[0] in gene_para_dict:
                    set_name = gene_para_dict[line[0]]
                    para_dict[set_name].AddGene(line[2], line[3])
                    gene_para_dict[line[2]] = set_name
    return para_dict, gene_para_dict

def FindOppoSet(id1_list, homo_dict, gene_para2, setToUnion2):
    oppo_set = set()
    assign_set = set()
    assign_union = set()
    for id1 in id1_list:
        # id1: gene in species1
        if id1 in homo_dict:
            # id1 has homologous genes
            oppo_idlist = homo_dict[id1].id2
            for id2 in oppo_idlist:
                if id2 in gene_para2:
                    # id2's paralogous genes
                    set2 = gene_para2[id2]
                    oppo_set.add(set2)
                    if set2 in setToUnion2:
                        assign_set.add(set2)
                        assign_union.add(setToUnion2[set2])
    return oppo_set, assign_set, assign_union

def UpdateSet2Union(assign_union, Union2Set, Set2Union):
    '''
    Redirect sets in SetToUnion.
    '''
    keep_union = assign_union[0]
    for un in assign_union[1:]:
        id_list = Union2Set[un]
        for id1 in id_list:
            Set2Union[id1] = keep_union
    return

def UpdateUnion2Set(assign_union, Union2Set):
    '''
    Combine unions, delete redundant unions.
    '''
    keep_union = assign_union[0]
    for un in assign_union[1:]:
        Union2Set[keep_union] += Union2Set[un]
        Union2Set.pop(un, None)
    return


def overlap(bed,tuple1):
    if bed.chr != tuple1[0]:
        return False
    if (bed.stop > tuple1[1]) and (bed.start < tuple1[2]):
        return True
    else:
        return False

def overlap_tuple(tuple1, tuple2):
    if tuple1[0] != tuple2[0]:
        return False
    if (tuple1[2] > tuple2[1]) and (tuple1[1] < tuple2[2]):
        return True
    return False


def CombineSet(set_list, para_dict):
    id_list = []
    name_list = []
    for set1 in set_list:
        id_list += list(para_dict[set1].id)
        name_list += list(para_dict[set1].name)
    id_list_comb = Unique(id_list)
    name_list_comb = Unique2(id_list, name_list)
    return id_list_comb, name_list_comb


def Main():
    args = ParseArg()
    cluster_n = 1

    # Paralogue files. Find paralogous gene cluster in each species.
    Para1, Gene_para1 = ParseParalog(args.paralog[0])
    Para2, Gene_para2 = ParseParalog(args.paralog[1])

    # Homologue files
    Ortho1 = ParseOrtho(args.ortho[0])
 #   Ortho2 = ParseOrtho(args.ortho[1])

    # Homologous clusters
    i = 1
    SetToUnion1 = {}
    SetToUnion2 = {}
    UnionToSet1 = {}
    UnionToSet2 = {}
    for set1 in Para1:
        #set1: set name in species1
        # oppo_set_list: a list of set names in species2. assign_set: a list of set names in species2 that are assigned to unions.
        # assign_union: a list of union names.
        oppo_set_list, assign_set, assign_union = FindOppoSet(Para1[set1].id, Ortho1, Gene_para2, SetToUnion2)
        oppo_set_list = list(oppo_set_list)
        assign_set = list(assign_set)
        assign_union = list(assign_union)

        if len(oppo_set_list) == 0:
            continue
        if len(assign_set) == 0:
            # the corresponding sets haven't been assigned to any unions
            union_name = "Union" + str(i)
            SetToUnion1[set1] = union_name
            UnionToSet1[union_name] = [set1]
            UnionToSet2[union_name] = []
            for set2 in oppo_set_list:
                SetToUnion2[set2] = union_name
                UnionToSet2[union_name].append(set2)
        elif len(assign_set) == 1:
            # one of the corresponding sets has been assigned to a union.
            print Para1[assign_set[0]].id
            union_name = SetToUnion2[assign_set[0]]
            SetToUnion1[set1] = union_name
            UnionToSet1[union_name].append(set1)
            for set2 in oppo_set_list:
                if set2 not in assign_set:
                    SetToUnion2[set2] = union_name
                    UnionToSet2[union_name].append(set2)
        else:
            # more than one corresponding sets have been assigned to different unions. Need to collapse.
            print Para1[assign_set[0]].id
            UpdateSet2Union(assign_union, UnionToSet2, SetToUnion2)
            UpdateUnion2Set(assign_union, UnionToSet2)

        i += 1
        if i % 1000 == 0:
            print >> sys.stderr, "Processed clusters:" + str(i)

    with open(args.output, "w") as fout:
        j = 1
        for un in UnionToSet1:
            id_set1, name_set1 = CombineSet(UnionToSet1[un], Para1)
            id_set2, name_set2 = CombineSet(UnionToSet2[un], Para2)
            if len(id_set1) > 1 or len(id_set2) > 1:
                for ID1, name1 in zip(id_set1, name_set1):
                    for ID2, name2 in zip(id_set2, name_set2):
                        print >>fout, "\t".join([ID1, name1, ID2, name2, "Cluster_" + str(j)])
            j += 1

Main()









