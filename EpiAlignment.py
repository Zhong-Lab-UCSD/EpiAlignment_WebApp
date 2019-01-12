# This is the version for EpiAlignment Web Service

import argparse
import sys
import copy
import math
from math import log, exp
from time import time
from multiprocessing import *
from functools import partial


class HomoRegion:
    '''
    class of homologous regions.
    '''

    def __init__(self):
        self.S1 = []
        self.S2 = []
        self.L = 0
        self.name = ""
        self.averagedL = 0
        self.loc2 = -1
        self.loc1 = -1
        self.start_point = 0
        self.prob = []
        self.S1_path = ""
        self.S2_path = ""


def ParseArg():
    p = argparse.ArgumentParser(description="EpiAlignment. A semi-global alignment algorithm for chromosomal similarity search.")
    p.add_argument("Input", type=str, help="Input file name.")
    p.add_argument(
        "-e",
        "--equil_file",
        type=str,
        help="The parameter file containing intial guesses for s, mu, k, equilibrium probabilities and weights.")
    p.add_argument("-p", "--process_num", type=int, default=1, help="Number of processes to be used. Default: 1.")
    p.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name. This file contains region name, alignment scores and target position for each region pair.")
    p.add_argument(
        "-O",
        "--out_allvec",
        type=str,
        help="Output file name. This file contains the last rows of the alignment matrices (alignment scores across the search regions). Only available when --all_prob is specified.")
    p.add_argument(
        "-r",
        "--align_path",
        type=str,
        help="Alignment path file name. The alignment path will be output if specified. WARNING: The reconstruction of alignment paths has excessive memory demand. Use only when input sequence number is small. ")
    if len(sys.argv) == 1:
        print >>sys.stderr, p.print_help()
        sys.exit(0)
    return p.parse_args()


def ReadInput(fin_name):
    '''
    Read the input file.
    fin_name: input file name, which is the the input parameter Input.
    return: a list of input region pairs. Each element is a HomoRegion object.
    '''
    Slist = []
    s1_count = 0
    s2_count = 0
    flag = 0
    with open(fin_name, "r") as fin:
        S = []
        s1_maxlen = 0
        s2_maxlen = 0
        s1_avelen = 0
        s2_avelen = 0
        line = fin.readline().strip()
        if "@" not in line:
            raise Exception(301, "The first line of the input file does not start with @.")
        flag = 1
        while True:
            line = fin.readline().strip()
            if len(line) == 0:
                if s1_count > s2_count:
                    Sobj.S2 = S
                    Slist.append(Sobj)
                    s2_count += 1
                    if len(Sobj.S1) > s1_maxlen:
                        s1_maxlen = len(Sobj.S1)
                    if len(Sobj.S2) > s2_maxlen:
                        s2_maxlen = len(Sobj.S2)
                    s1_avelen += len(Sobj.S1)
                    s2_avelen += len(Sobj.S2)
                else:
                    raise Exception(302, "The number of sequences are different!")
                break
            if line == "+":
                i = 0
                flag = 0
                continue
            if "@" in line:
                if s1_count > s2_count:
                    Sobj.S2 = S
                    Slist.append(Sobj)
                    s2_count += 1
                    if len(Sobj.S1) > s1_maxlen:
                        s1_maxlen = len(Sobj.S1)
                    if len(Sobj.S2) > s2_maxlen:
                        s2_maxlen = len(Sobj.S2)
                    s1_avelen += len(Sobj.S1)
                    s2_avelen += len(Sobj.S2)
                else:
                    Sobj = HomoRegion()
                    Sobj.name = line[1:]
                    Sobj.S1 = S
                    s1_count += 1
                flag = 1
                S = []
                continue
            if flag == 1:
                line = line.upper()
                S += [(x, "") for x in line]
            else:
                S = S[0:i] + [(a[0], a[1] + b) for a, b in zip(S[i:(i + len(line))], line)] + S[(i + len(line)):]
                i += len(line)
    return Slist, s1_maxlen, s2_maxlen, s1_avelen / len(Slist), s2_avelen / len(Slist)


def ReadParameters(f_name):
    '''
    Read parameters. Build the dictionaries of equilibrium probabilities on the linear and log scales.
    f_name: parameter file name, which is specified by --equil_file.
    Sample return:
    parameter vector x: [0.1, 0.01, 0.1]
    weights: [1.0, 0.0]
    equil_dic: {'A': 0.25, 1: [0.9, 0.1], 'C': 0.25, 'T': 0.25, 'G': 0.25}
    log_equil_dict: {'A': -1.386, 1: [-0.105, -2.303], 'C': -1.386, 'T': -1.386, 'G': -1.386}
    '''
    n_epi = 0
    x = []
    weights = []
    equil_dict = {}
    log_equil_dict = {}
    with open(f_name, "r") as fin:
        x.append(float(fin.readline().strip()))
        x.append(float(fin.readline().strip()))
        for line in fin:
            line = line.strip().split("\t")
            if len(line) > 1:
                if n_epi == 0:
                    for item in line:
                        equil_dict[item.split(":")[0]] = float(item.split(":")[1])
                        log_equil_dict[item.split(":")[0]] = log(float(item.split(":")[1]))
                    n_epi += 1
                elif n_epi > 0:
                    if ":" in line[0]:
                        equil_dict[n_epi] = [0.0, 0.0]
                        log_equil_dict[n_epi] = [0.0, 0.0]
                        for item in line:
                            equil_dict[n_epi][int(item.split(":")[0])] = float(item.split(":")[1])
                            log_equil_dict[n_epi][int(item.split(":")[0])] = log(float(item.split(":")[1]))
                        n_epi += 1
                    else:
                        weights = [float(n) for n in line]
            else:
                x.append(float(line[0]))
    return x, weights, equil_dict, log_equil_dict

# def Log_sum(lnA, lnB, lnC, lnD):
#   '''
#   tmp=[logA, logB, logC, logD]
#   return log(A+B+C-D)
#   '''

#   tm = max(lnA, lnB, lnC, lnD)
#   sm = exp(lnA - tm) + exp(lnB - tm) + exp(lnC - tm) - exp(lnD - tm)

#   try:
#     return tm + log(sm)
#   except:
#     return float("-Inf")

# def Log_sum3(lnA, lnB, lnC):
#   '''
#   only for path_len_mat
#   tmp=[logA, logB, logC]
#   return log(A+B+C)
#   '''

#   tm = max(lnA, lnB, lnC)
#   sm = exp(lnA - tm) + exp(lnB - tm) + exp(lnC - tm)

#   try:
#     return tm + log(sm)
#   except:
    return float("-Inf")


# def Path_matrix(mlen1, mlen2, mu, lamb, beta, link_p, log_link_p):
#     '''
#     Pre-compute the path matrix for normalization.
#     mlen1, mlen2: maximal lengths of ancestral and descendent regions.
#     mu, lamb, beta, link_p: parameters
#     return: a mlen1 * mlen2 matrix.
#     '''
#     if mlen1 > mlen2:
#         tmp = mlen2
#         mlen2 = mlen1
#         mlen1 = tmp
#     if mlen1 > mlen2:
#         mlen2 = mlen1
#     else:
#         mlen1 = mlen2
#     path_mat3 = [[Na] * (mlen2 + 1) for i in xrange(mlen1 + 1)]
#     path_mat2 = [[Na] * (mlen2 + 1) for i in xrange(mlen1 + 1)]

#     link_p12 = max(link_p[1], link_p[2])
#     log_link_p12 = log(link_p12)
#     log_lamb_mu = log(lamb / mu)
#     lamb_beta = lamb * beta
#     log_lamb_beta = log(lamb_beta)

    
#     path_mat3[0][0] = log_link_p[3] + log(Gamma(0,lamb,mu))
#     path_mat3[1][0] = log_link_p[3] + log_link_p[0] + log(Gamma(0,lamb,mu)) + log_lamb_mu

#     path_mat3[0][1] = log_link_p[3] + log_lamb_beta + log(Gamma(0,lamb,mu))
#     path_mat2[0][1] = log_link_p[3] + log_lamb_beta + log(Gamma(0,lamb,mu))

#     # First row
#     for i in xrange(2, mlen2 + 1):
#         path_mat3[0][i] = path_mat3[0][i-1] + log_lamb_beta
#         path_mat2[0][i] = path_mat2[0][i-1] + log_lamb_beta
#     # First column
#     for i in xrange(2, mlen1 + 1):
#         path_mat3[i][0] = path_mat3[i-1][0] + log_link_p[0] + log_lamb_mu
#     for i in xrange(1, (mlen1 + 1)):
#         for j in xrange(1, (mlen2 + 1)):

#             ent0 = log_lamb_mu + log_link_p[0] + path_mat3[i - 1][j]
#             ent1 = log_link_p12 + log_lamb_mu + path_mat3[i - 1][j - 1]
#             ent2 = log_lamb_beta + path_mat2[i][j - 1]

#             max_v3, max_v2, tup_ind1, tup_ind2 = Maximum(ent0, ent1, ent2, j)
#             #print "max"
#             #print max_v3, max_v2
#             #print path_mat3
#             #print path_mat2
#             #print "\n"

#             path_mat3[i][j] = max_v3
#             path_mat2[i][j] = max_v2


# #    for line in path_mat3:
# #        print " ".join([str(f) for f in line])
# #    print "\n"
#     return path_mat3



# def Mod_equilibrium(e_dict, log_e_dict, weights):
#     '''
#     Scale the equilibrium probabilities for DNA bases by sequence weight.
#     e_dict, log_e_dict: dictionaries of equilibrium probabilities on the linear and log scales. equil_dict and log_equil_dict returned by ReadParameters.
#     weights: the weights vector returned by ReadParameters.
#     return: None. The function will modify the dictionaries directly.
#     '''
#     for b in "A", "C", "G", "T":
#         e_dict[b] = e_dict[b]**weights[0]
#         log_e_dict[b] = log_e_dict[b] * weights[0]


def Epi_equilibrium(n_epi, equil_dict, log_equil_dict, weights):
    '''
    Products of the equilibrium probabilities of epigenetic marks.
    The function will iterate all possible combinations of '1's and '0's (k) and compute the equilibrium probability of observing the combination.
    For example, if there are two epi-marks, all possible k's will be '00', '01', '10', '11'.
    The equilibrium probabilities are scaled by epi-weights.
    n_epi: number of epigenomic marks.
    weights: the weights vector.
    return: two dictionaries, in which the keys are combination of '1's and '0's and values are equilibrium probabilities.
    '''
    S_epi = {}
    log_S_epi = {}
    for i in xrange(pow(2, n_epi)):
        k = bin(i)[2:].zfill(n_epi)
        v = 1.0
        lv = 0.0
        for j in xrange(0, n_epi):
            v = v * (equil_dict[j + 1][int(k[j])]) ** weights[j + 1]
            lv += weights[j + 1] * log_equil_dict[j + 1][int(k[j])]
        S_epi[k] = v
        log_S_epi[k] = lv

    return S_epi, log_S_epi


def Equilibrium_matrix(log_equil_dict, log_S_epi, weights):
    '''
    Products of the equilibrium probabilities of bases and epigenetic marks.
    The base equilibrium probabilities are scaled by sequence-weights.
    return: a dictionary. Keys: base - epi-state. Values: scaled equilibrium probabilities.
    '''
    log_equil_mat = {}
    for base in "A", "C", "G", "T":
        log_equil_mat[base] = {}
        for epi in log_S_epi:
            log_equil_mat[base][epi] = log_equil_dict[base] * weights[0] + log_S_epi[epi]
    return log_equil_mat


def Link_prob(prime, n, b, lam, mu):
    """
    Compute p', p'', p''' defined in the TKF DNA evolutionary model.
    prime: the prime of p (0, 1 or 2)
    n: the subscript of p.
    b, lam, mu: the value of beta, lambda and mu. See the function Manhattan for the definition of beta.
    return: p'_0, p_1, p'_1 or p''_1
    """
    if prime == 1 and n == 0:
        return mu * b
    elif prime == 0 and n == 1:
        return exp(-1 * mu) * (1 - lam * b)
    elif prime == 1 and n == 1:
        return (1 - exp(-1 * mu) - mu * b) * (1 - lam * b)
    elif prime == 2 and n == 1:
        return (1 - lam * b)


def Transition_f(s, base1, base2, e_dict):
    '''
    The transition probability function f between DNA bases, defined in the TKF DNA evolutionary model.
    s: DNA base substitution rate.
    base1: a DNA base in the query sequence (A, C, G or T).
    base2: a DNA base in the search sequence (A, C, G or T).
    e_dict: the dictionary of equilibrium probability.
    return: transition probability between base 1 and 2, which is the value of function f.
    '''
    e = exp(-s)
    if base1 == base2:
        return e + e_dict[base2] * (1 - e)
    else:
        return e_dict[base2] * (1 - e)


def Transition_g(i, e1, e2, k, equil_dict):
    '''
    The transition probability function g between epigenomic states.
    i: the index of epi mark (the i-th epi mark).
    e1: an epi state in the query region ('1' or '0').
    e2: an epi state in the search sequence ('1' or '0').
    k: epigenomic state changing rate.
    return: transition probability between epi states 1 and 2, which is the value of function g.
    '''
    e = exp(-k)
    if e1 == e2:
        return e + equil_dict[i][int(e2)] * (1 - e)
    else:
        return equil_dict[i][int(e2)] * (1 - e)


def Trans_matrix(n_epi, x, e_dict, weights):
    '''
    Construct the dictionary of transition probabilities.

    n_epi: number of epi marks.
    x: the parameter vector with s, mu and kappas.
    e_dict: the dictionary of equilibrium probabilities.
    weights: the weights vector.

    return: a dictionary. For the four bases, the keys are the four bases. For each base, the value is a dict,
    in which the keys are the four bases, and the values are transition probabilities between bases on log scales.
    For epigenomic marks, the keys are the indices of epi marks. For each index, the value is a dict,
    in which the keys are '1' and '0', and the values are dicts with '1' and '0' as keys and transition probabilities between epi states on log scales as values.
    Keys: the i-th epigenetic mark.
    Sample format: {i:{'0':{'0':a,'1':b},'1':{'0':c,'1':d}}, 'A': {'A':aa, 'C':ac, 'G':ag, 'T':at}}
    '''
    log_trans_dic = {}
    base = ["A", "C", "G", "T"]
    for b in base:
        log_trans_dic[b] = {}
        for b1 in base:
            log_trans_dic[b][b1] = log(Transition_f(x[0], b, b1, e_dict)) * weights[0]
    for i in xrange(1, n_epi + 1):
        log_trans_dic[i] = {}
        for e1 in ['0', '1']:
            log_trans_dic[i][e1] = {}
            for e2 in ['0', '1']:
                log_trans_dic[i][e1][e2] = log(Transition_g(i, e1, e2, x[1 + i], e_dict)) * weights[i]
    return log_trans_dic

def Combine_epi_trans(log_trans_dic, log_S_epi):
    '''
    Combine the epi trans matrix.
    '''
    log_trans_prod = {}
    for epi1 in log_S_epi:
        log_trans_prod[epi1] = {}
        for epi2 in log_S_epi:
            g = 0
            for i in xrange(1, len(epi1) + 1):
                g += log_trans_dic[i][epi1[i - 1]][epi2[i - 1]]
            log_trans_prod[epi1][epi2] = g
    return log_trans_prod



# def Transition_g_sum(tuple1, tuple2, log_trans_dic):
#     '''
#     Compute the transition probabilities between two sets of epigenomic states.
#     tuple1, 2: the epigenomic states at a specific position (for example, the m-th position in S1 and the n-th position in S2).
#     log_trans_dic: the dictionary containing transition probabilities on log scale, constructed by Trans_matrix.
#     '''
#     g = 0.0
#     if hypN == 1 or tuple1[0] == tuple2[0]:
#         for i in xrange(1, len(tuple1)):
#             g = g + log_trans_dic[i][tuple1[i]][tuple2[i]]
#         return g
#     else:
#         return log_S_epi["".join(tuple2[1:])]


def Backtrack_mat(ind, ind2, bt, i, j, sp):
    '''
    Update the backtrack matrix and record the alignment path.
    ind, ind2: indices of maxima.
    bt: the backtrack matrix.
    i, j: current row and column indices.
    sp: start point. Either 0 or 1.
    '''
    ind += sp
    if ind == 0:
        bt[i][j] = "u"
    elif ind == 1:
        if ind2 == 0:
            bt[i][j] = "d"
        elif ind2 == 1:
            bt[i][j] = "z"
    elif ind == 2:
        bt[i][j] = "l"


def Recons_path(S, bt, j_start, i_start):
    '''
    Reconstruct the optimal alignment path from the backtrack matrix.
    S: a HomoRegion object.
    bt: the backtrack matrix.
    j_start: the index at which the maximal alignment score is found.
    '''
    if len(S.S1) <= len(S.S2):
        S1 = S.S1
        S2 = S.S2
    else:
        S1 = S.S2
        S2 = S.S1

    n_epi = len(S.S1[0]) - 1
    S_match = ""
    S1_epi_path = {}
    S2_epi_path = {}
    for k in xrange(1, n_epi + 1):
        S1_epi_path[k] = ""
        S2_epi_path[k] = ""

    i = i_start  # S1
    j = j_start  # S2

    S1_align = ""
    S2_align = ""
    while i != 0 and j != 0:
        if bt[i][j] == "d":
            S1_align += S1[i - 1][0]
            S2_align += S2[j - 1][0]
            for k in xrange(1, n_epi + 1):
                S1_epi_path[k] += S1[i - 1][k]
                S2_epi_path[k] += S2[j - 1][k]
            if S1[i - 1][0] == S2[j - 1][0]:
                S_match += "|"
            else:
                S_match += " "
            i = i - 1
            j = j - 1
        elif bt[i][j] == "l":
            S1_align += "-"
            S2_align += S2[j - 1][0]
            for k in xrange(1, n_epi + 1):
                S1_epi_path[k] += "-"
                S2_epi_path[k] += S2[j - 1][k]
            S_match += " "
            j = j - 1
        elif bt[i][j] == "u":
            S1_align += S1[i - 1][0]
            S2_align += "-"
            for k in xrange(1, n_epi + 1):
                S1_epi_path[k] += S1[i - 1][k]
                S2_epi_path[k] += "-"
            S_match += " "
            i = i - 1
        elif bt[i][j] == "z":
            S1_align += (S1[i - 1][0] + "-")
            S2_align += ("-" + S2[j - 1][0])
            for k in xrange(1, n_epi + 1):
                S1_epi_path[k] += (S1[i - 1][k] + "-")
                S2_epi_path[k] += ("-" + S2[j - 1][k])
            S_match += " "
            i = i - 1
            j = j - 1

    for k in xrange(1, n_epi + 1):
        S1_epi_path[k] = S1_epi_path[k][::-1]
        S2_epi_path[k] = S2_epi_path[k][::-1]

    if len(S.S1) < len(S.S2):
        S.S1_path = S1_align[::-1].lstrip("-")
        S.S2_path = S2_align[len(S.S1_path) - 1::-1]
        S.S_match = S_match[::-1]
        S.S1_epi_path = S1_epi_path
        S.S2_epi_path = S2_epi_path
    else:
        S.S2_path = S1_align[::-1].lstrip("-")
        S.S1_path = S2_align[len(S.S1_path) - 1::-1]
        S.S_match = S_match[::-1]
        S.S2_epi_path = S1_epi_path
        S.S1_epi_path = S2_epi_path

    return


def Print_path(S0, fout2):
    '''
    Print out the alignment path when "--align_path" is specified.
    S0: a HomoRegion.
    fout2: the output file where to print the alignment path.
    '''
    n_epi = len(S0.S1[0]) - 1
    j = 0
    while j + 100 < len(S0.S1_path):
        print >>fout2, "@Sequence name: " + S0.name

        for i in range(1, n_epi + 1)[::-1]:
            print >>fout2, S0.S1_epi_path[i][j:j + 100]

        print >>fout2, S0.S1_path[j:j + 100]
        print >>fout2, S0.S_match[j:j + 100]
        print >>fout2, S0.S2_path[j:j + 100]

        for i in xrange(1, n_epi + 1):
            print >>fout2, S0.S2_epi_path[i][j:j + 100]
        print >>fout2, ""
        j += 100

    print >>fout2, "@Sequence name: " + S0.name

    for i in range(1, n_epi + 1)[::-1]:
        print >>fout2, S0.S1_epi_path[i][j:j + 100]

    print >>fout2, S0.S1_path[j:j + 100]
    print >>fout2, S0.S_match[j:j + 100]
    print >>fout2, S0.S2_path[j:j + 100]

    for i in xrange(1, n_epi + 1):
        print >>fout2, S0.S2_epi_path[i][j:j + 100]
    print >>fout2, ""


def Maximum(m0, m1, m2, j):
    '''
    Return the maximal value in m0, m1 and m2 and the maximal value in m1, m2.
    '''

    if m1 >= m2:
        # manh1>manh2
        if m0 >= m1:
            # manh0>manh1
            return m0, m1, 0, j
        else:
            # manh1 is the largest
            return m1, m1, 0, (j - 1)
    else:
        # manh2>manh1
        if m0 >= m2:
            #manh0 > manh2
            return m0, m2, 0, j
        else:
            # manh2 is the largest
            return m2, m2, 1, (j - 1)



def Gamma(n, lam, mu):
  '''
  Equilibrium probability of a sequence of length n
  '''
  return (1 - lam / mu) * math.pow(lam / mu,n)



def Path_norm_delta(start_point, prev_point, current_point):
    '''
    Find the normalization delta.
    start_point, prev_point, current_point: tuples. 
    return: delta. path_nmat[prev_point - start_point] - path_nmat[current_point - start_point]
    '''
    prev_delta = max(prev_point[0] - start_point[0], prev_point[1] - start_point[1])
    current_delta = max(current_point[0] - start_point[0], current_point[1] - start_point[1])
    if prev_delta == current_delta:
        return 0
    else:
        return (- diag_norm)

def Manhattan(S, X):
    '''
    Initialize and fill the matrices in dynamic programming for alignment score computation.
    S: a HomoRegion object.
    X: the parameter vector containing s, mu and kappas.
    log_trans_dic: the dictionary containing transition probabilities on log scale, constructed by Trans_matrix.
    Note that the argument to be distributed to different processes should be the first one.
    return: the updated S.
    '''
    if len(S.S1) <= len(S.S2):
        S1 = S.S1
        S2 = S.S2
    else:
        S1 = S.S2
        S2 = S.S1
    m = len(S1)
    n = len(S2)
    s = X[0]
    kappa = X[2:]

    # Initialization
    if align_path:
        # backtrack matrix: upper, left, diagonal
        backtrack0 = [[0] * (n + 1) for i in xrange(m + 1)]
        backtrack0[0] = ["l"] * (n + 1)
        for b in range(1, m + 1):
            backtrack0[b][0] = "u"

    # maximum of three values (manh0, manh1, manh2), maximum of two values (manh1, manh2)
    manh3 = [[Na] * (n + 1) for i in xrange(2)]
    manh2 = [[Na] * (n + 1) for i in xrange(2)]
    # Start point. The first row and column: all start from the position itself.
    manh3_st = [[Na] * (n + 1) for i in xrange(2)]
    # last column
    last_col = [Na] * (m + 1)
    # last colum start point
    last_col_st = [Na] * (m + 1)

    # 0,0
    init0 = 0
    # init0 = log_link_p[3] + log(Gamma(0,lamb,mu))

    for i in xrange(0, n + 1):
        manh3[0][i] = init0
        manh2[0][i] = init0
        manh3_st[0][i] = (0, i)

    # 1,0
    manh3[1][0] = init0
    manh3_st[1][0] = (1, 0)

    ent0_comp = log_lamb_mu + log_link_p[0]

    # Dynamic programming starts
    for i in xrange(1, (m + 1)):
        for j in xrange(1, (n + 1)):
            tmp0 = Log_transition_dic[S1[i - 1][0]][S2[j - 1][0]] + log_link_p[1] + Log_trans_prod[S1[i - 1][1]][S2[j - 1][1]]
            tmp1 = log_link_p[2] + log_equil_mat[S2[j - 1][0]][S2[j - 1][1]]

            max_t = max(tmp0, tmp1)

            ent0 = ent0_comp + manh3[0][j] - half_diag_norm

            ent2 = log_lamb_beta + manh2[1][j - 1] - half_diag_norm

            ent1 = log_lamb_mu + max_t + manh3[0][j - 1] - log_equil_mat[S2[j - 1][0]][S2[j - 1][1]] - diag_norm

            max_v3, max_v2, tup_ind1, tup_ind2 = Maximum(ent0, ent1, ent2, j)

            manh3[1][j] = max_v3
            manh2[1][j] = max_v2
            manh3_st[1][j] = manh3_st[tup_ind1][tup_ind2]

            if align_path:
                max_it = tmp.index(max_t)
                index_v3 = [ent0, ent1, ent2].index(max_v3)
                Backtrack_mat(index_v3, max_it, backtrack0, i, j, 0)


        last_col[i] = manh3[1][-1]
        last_col_st[i] = manh3_st[1][-1]

        manh3[0] = manh3[1]
        manh3[1] = [Na] * (n + 1)

        manh2[0] = manh2[1]
        manh2[1] = [Na] * (n + 1)

        manh3_st[0] = manh3_st[1]
        manh3_st[1] = [Na] * (n + 1)
        #print manh3_st[0]

        if i < m:
            manh3[1][0] = init0

            manh3_st[1][0] = (i + 1, 0)


  #  print norm_max_list
  #  print plen_flist
  #  print last_col
    last_row_max = max(manh3[0][1:])
    last_col_max = max(last_col[1:])


    S.L = max(last_row_max, last_col_max)
    if last_row_max >= last_col_max:
        i_start = m
        j_start = manh3[0].index(last_row_max)
        S.start_point = manh3_st[0][j_start]
    else:
        i_start = last_col.index(last_col_max)
        j_start = n
        S.start_point = last_col_st[i_start]

    #for line in backtrack0:
    #    print >> sys.stderr, " ".join(line)
    S.loc2 = j_start
    S.loc1 = i_start
    S.averagedL = sum(manh3[0][1:] + last_col[1:]) / float(m + n)
    if all_prob:
        S.prob = manh3[0][1:] + last_col[1:]
    if align_path:
        Recons_path(S, backtrack0, j_start, i_start)
    return S


def Manhattan_obj(x, S):
    '''
    Distribute region pairs to difference processes for parallel computing.
    x: the parameter vector.
    S: the list of region pairs. Each element is a HomoRegion object.
    e_dict: the dictionary of equilibrium probabilities.
    weights: the weights vector.
    return: None. The function will update the S list directly.
    '''
    try:
        p = Pool(p_num)
        mp_queue = p.map(partial(Manhattan, X=x), S)
        S[:] = mp_queue
        p.close()
        p.join()
    except Exception as err:
        p.terminate()
        print >> sys.stderr, err.args[1]
        sys.exit(err.args[0])
    # Manhattan(S[0], x)
    return


def Main():
    args = ParseArg()
    S, maxlen1, maxlen2, ave1, ave2 = ReadInput(args.Input)

    global hypN
    hypN = 1

    global p_num
    p_num = args.process_num

    global all_prob
    all_prob = args.out_allvec

    global align_path
    align_path = args.align_path

    x, weights, equil_dict, log_equil_dict = ReadParameters(args.equil_file)

    global Na
    Na = float('-Inf')

    global mu, lamb, beta, log_link_p
    mu = x[1]
    lamb = mu * (ave1 + ave2) / (ave1 + ave2 + 2)
    beta = (1 - exp(lamb - mu)) / (mu - lamb * exp(lamb - mu))
    if (1 - exp(-mu) - mu * beta) < 0:
        raise Exception(206, "Invalid parameters.")
    link_p = [Link_prob(1, 0, beta, lamb, mu), Link_prob(0, 1, beta, lamb, mu), Link_prob(1, 1, beta, lamb, mu), Link_prob(2, 1, beta, lamb, mu)]
    log_link_p = [log(p) for p in link_p]

    global diag_norm, half_diag_norm, log_lamb_mu, log_lamb_beta
    log_lamb_mu = log(lamb / mu)
    log_lamb_beta = log(lamb * beta)
    diag_norm = max(log_link_p[1], log_link_p[2]) + log_lamb_mu
    half_diag_norm = 0.5 * diag_norm

    # Equilibrium probabilities
    global log_equil_mat
    S_epi, log_S_epi = Epi_equilibrium(len(S[0].S1[0]) - 1, equil_dict, log_equil_dict, weights)
    log_equil_mat = Equilibrium_matrix(log_equil_dict, log_S_epi, weights)

    # Transition_matrix
    global Log_transition_dic, Log_trans_prod
    Log_transition_dic = Trans_matrix(len(S[0].S1[0]) - 1, x, equil_dict, weights)
    Log_trans_prod = Combine_epi_trans(Log_transition_dic, log_S_epi)


    # t0 = time()
    Manhattan_obj(x, S)
    # t1 = time()
    if all_prob:
        with open(args.output, "w") as fout, open(args.out_allvec, "w") as fout2:
            for pair in S:
                print >>fout, "\t".join([pair.name, str(pair.L), str(pair.averagedL), str(pair.start_point[0]), str(pair.loc1), str(pair.start_point[1]), str(pair.loc2)])
                print >>fout2, ",".join([str(f) for f in [pair.name] + pair.prob])
    else:
        with open(args.output, "w") as fout:
            for pair in S:
                print >>fout, "\t".join([pair.name, str(pair.L), str(pair.averagedL), str(pair.start_point[0]), str(pair.loc1), str(pair.start_point[1]), str(pair.loc2)])

    if align_path:
        with open(align_path, "w") as fout2:
            for pair in S:
                Print_path(pair, fout2)

    # print >>sys.stderr, "Time:%f" % ((t1 - t0) / 60)


Main()
