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
        self.loc = -1
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
        print(p.print_help(), file=sys.stderr)
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
        line = fin.readline().strip()
        if "@" not in line:
            print("The input file is not FASTQ", file=sys.stderr)
            exit()
        flag = 1
        while True:
            line = fin.readline().strip()
            if len(line) == 0:
                if s1_count > s2_count:
                    Sobj.S2 = S
                    Slist.append(Sobj)
                    s2_count += 1
                else:
                    print("The number of sequences are different!", file=sys.stderr)
                    exit()
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
                S += [(x,) for x in line]
            else:
                S = S[0:i] + [a + (b,) for a, b in zip(S[i:(i + len(line))], line)] + S[(i + len(line)):]
                i += len(line)
    return Slist


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


def Mod_equilibrium(e_dict, log_e_dict, weights):
    '''
    Scale the equilibrium probabilities for DNA bases by sequence weight.
    e_dict, log_e_dict: dictionaries of equilibrium probabilities on the linear and log scales. equil_dict and log_equil_dict returned by ReadParameters.
    weights: the weights vector returned by ReadParameters.
    return: None. The function will modify the dictionaries directly.
    '''
    for b in "A", "C", "G", "T":
        e_dict[b] = e_dict[b]**weights[0]
        log_e_dict[b] = log_e_dict[b] * weights[0]


def Epi_equilibrium(e_dict, log_e_dict, n_epi, weights):
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
    for i in range(pow(2, n_epi)):
        k = bin(i)[2:].zfill(n_epi)
        v = 1.0
        lv = 0.0
        for j in range(0, n_epi):
            v = v * (e_dict[j + 1][int(k[j])]) ** weights[j + 1]
            lv += weights[j + 1] * log_e_dict[j + 1][int(k[j])]
        S_epi[k] = v
        log_S_epi[k] = lv

    return S_epi, log_S_epi


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
        return math.exp(-1 * mu) * (1 - lam * b)
    elif prime == 1 and n == 1:
        return (1 - math.exp(-1 * mu) - mu * b) * (1 - lam * b)
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
    e = math.exp(-s)
    if base1 == base2:
        return e + e_dict[base2] * (1 - e)
    else:
        return e_dict[base2] * (1 - e)


def Transition_g(i, e1, e2, k, e_dict):
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
        return e + e_dict[i][int(e2)] * (1 - e)
    else:
        return e_dict[i][int(e2)] * (1 - e)


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
    for i in range(1, n_epi + 1):
        log_trans_dic[i] = {}
        for e1 in ['0', '1']:
            log_trans_dic[i][e1] = {}
            for e2 in ['0', '1']:
                log_trans_dic[i][e1][e2] = log(Transition_g(i, e1, e2, x[1 + i], e_dict)) * weights[i]
    return log_trans_dic


def Transition_g_sum(tuple1, tuple2, log_trans_dic, log_S_epi, hypN):
    '''
    Compute the transition probabilities between two sets of epigenomic states.
    tuple1, 2: the epigenomic states at a specific position (for example, the m-th position in S1 and the n-th position in S2).
    log_trans_dic: the dictionary containing transition probabilities on log scale, constructed by Trans_matrix.
    '''
    g = 0.0
    if hypN == 1 or tuple1[0] == tuple2[0]:
        for i in range(1, len(tuple1)):
            g = g + log_trans_dic[i][tuple1[i]][tuple2[i]]
        return g
    else:
        return log_S_epi["".join(tuple2[1:])]


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


def Recons_path(S, bt, j_start):
    '''
    Reconstruct the optimal alignment path from the backtrack matrix.
    S: a HomoRegion object.
    bt: the backtrack matrix.
    j_start: the index at which the maximal alignment score is found.
    '''
    if len(S.S1) < len(S.S2):
        S1 = S.S1
        S2 = S.S2
    else:
        S1 = S.S2
        S2 = S.S1

    n_epi = len(S.S1[0]) - 1
    S_match = ""
    S1_epi_path = {}
    S2_epi_path = {}
    for k in range(1, n_epi + 1):
        S1_epi_path[k] = ""
        S2_epi_path[k] = ""

    i = len(bt) - 1  # S1
    j = j_start  # S2

    S1_align = ""
    S2_align = ""
    while i != 0:
        if bt[i][j] == "d":
            S1_align += S1[i - 1][0]
            S2_align += S2[j - 1][0]
            for k in range(1, n_epi + 1):
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
            for k in range(1, n_epi + 1):
                S1_epi_path[k] += "-"
                S2_epi_path[k] += S2[j - 1][k]
            S_match += " "
            j = j - 1
        elif bt[i][j] == "u":
            S1_align += S1[i - 1][0]
            S2_align += "-"
            for k in range(1, n_epi + 1):
                S1_epi_path[k] += S1[i - 1][k]
                S2_epi_path[k] += "-"
            S_match += " "
            i = i - 1
        elif bt[i][j] == "z":
            S1_align += (S1[i - 1][0] + "-")
            S2_align += ("-" + S2[j - 1][0])
            for k in range(1, n_epi + 1):
                S1_epi_path[k] += (S1[i - 1][k] + "-")
                S2_epi_path[k] += ("-" + S2[j - 1][k])
            S_match += " "
            i = i - 1
            j = j - 1

    for k in range(1, n_epi + 1):
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
        print("@Sequence name: " + S0.name, file=fout2)

        for i in range(1, n_epi + 1)[::-1]:
            print(S0.S1_epi_path[i][j:j + 100], file=fout2)

        print(S0.S1_path[j:j + 100], file=fout2)
        print(S0.S_match[j:j + 100], file=fout2)
        print(S0.S2_path[j:j + 100], file=fout2)

        for i in range(1, n_epi + 1):
            print(S0.S2_epi_path[i][j:j + 100], file=fout2)
        print("", file=fout2)
        j += 100

    print("@Sequence name: " + S0.name, file=fout2)

    for i in range(1, n_epi + 1)[::-1]:
        print(S0.S1_epi_path[i][j:j + 100], file=fout2)

    print(S0.S1_path[j:j + 100], file=fout2)
    print(S0.S_match[j:j + 100], file=fout2)
    print(S0.S2_path[j:j + 100], file=fout2)

    for i in range(1, n_epi + 1):
        print(S0.S2_epi_path[i][j:j + 100], file=fout2)
    print("", file=fout2)


def Maxima(m_list, pl_list):
    '''
    m_list: manh0, manh1, manh2
    pl_list: path lengths
    The m_list will be normalized by the pl_list. The index of the maximum of the normalized list, mi3,
    as well as the index of the maximum of the last two elements, mi2 will be found.
    return: m_list[mi3], m_list[mi2], mi3, mi2+1
    '''
    plen_list = [(1.0 if f == 0 else f) for f in pl_list]
    norm_list = [f1 / f2 for f1, f2 in zip(m_list, plen_list)]

    if norm_list[1] >= norm_list[2]:
        # manh1>manh2
        if norm_list[0] >= norm_list[1]:
            # manh0>manh1
            return m_list[0], m_list[1], 0, 1
        else:
            # manh1 is the largest
            return m_list[1], m_list[1], 1, 1
    else:
        # manh2<manh1
        if norm_list[0] >= norm_list[2]:
            #manh0 > manh2
            return m_list[0], m_list[2], 0, 2
        else:
            # manh2 is the largest
            return m_list[2], m_list[2], 2, 2


def LocalPi(S1, log_trans_dic):
    '''
    Compute local equilibrium probabilities.
    S1: a genomic region, the S1 or S2 attribute of a HomoRegion object.
    log_trans_dic: the dictionary containing transition probabilities on log scale, constructed by Trans_matrix.
    return: a dictionary, in which the keys are combinations of '1's and '0's, and the values are local equilibrium probabilities.
    '''
    epi_n = len(S1[0][1:])
    log_SS_epi = {}
    phi_list = []
    m = len(S1) * 1.0

    for i in range(1, epi_n + 1):
        s = sum([1 for b in S1 if b[i] == '1'])
        phi_list.append([1 - s / m, s / m])

    for i in range(pow(2, epi_n)):
        k = bin(i)[2:].zfill(epi_n)
        v = 1.0
        for j in range(0, epi_n):
            v = v * phi_list[j][0] * exp(log_trans_dic[j + 1]['0'][k[j]]) + phi_list[j][1] * exp(log_trans_dic[j + 1]['1'][k[j]])
        log_SS_epi[k] = log(v)

    return log_SS_epi


def Manhattan(S, X, log_trans_dic, log_equil_dict, params):
    '''
    Initialize and fill the matrices in dynamic programming for alignment score computation.
    S: a HomoRegion object.
    X: the parameter vector containing s, mu and kappas.
    log_trans_dic: the dictionary containing transition probabilities on log scale, constructed by Trans_matrix.
    Note that the argument to be distributed to different processes should be the first one.
    return: the updated S.
    '''
    Na = float('-Inf')
    if len(S.S1) < len(S.S2):
        S1 = S.S1
        S2 = S.S2
    else:
        S1 = S.S2
        S2 = S.S1
    m = len(S1)
    n = len(S2)
    s = X[0]
    mu = X[1]
    kappa = X[2:]
    lamb = mu * (n + m) / (n + m + 2)

    # Initialization
    beta = (1 - math.exp(lamb - mu)) / (mu - lamb * math.exp(lamb - mu))
    link_p = [Link_prob(1, 0, beta, lamb, mu), Link_prob(0, 1, beta, lamb, mu), Link_prob(1, 1, beta, lamb, mu), Link_prob(2, 1, beta, lamb, mu)]
    if (1 - exp(-mu) - mu * beta) < 0:
        return float('Inf')

    log_link_p = [log(lp) for lp in link_p]

    log_lamb_mu = log(lamb / mu)
    log_lamb_beta = log(lamb * beta)

    log_S1_epi = LocalPi(S1, log_trans_dic)
    log_S2_epi = LocalPi(S2, log_trans_dic)
    log_S_epi = params['log_S_epi']

    if params['align_path']:
        # backtrack matrix: upper, left, diagonal
        backtrack0 = [[0] * (n + 1) for i in range(m + 1)]
        backtrack0[0] = ["l"] * (n + 1)
        for b in range(1, m + 1):
            backtrack0[b][0] = "u"

    # maximum of three values (manh0, manh1, manh2), maximum of two values (manh1, manh2)
    manh3 = [[Na] * (n + 1) for i in range(2)]
    manh2 = [[Na] * (n + 1) for i in range(2)]
    # path lengths
    manh3_plen = [[0] * (n + 1) for i in range(2)]
    manh2_plen = [[0] * (n + 1) for i in range(2)]

    # 0,0
    ent0 = 0

    for i in range(0, n + 1):
        manh3[0][i] = ent0
        manh2[0][i] = ent0

    # 1,0
    manh3[1][0] = ent0 + log_lamb_mu + log_equil_dict[S1[0][0]] + log_link_p[0] + log_S2_epi["".join(S1[0][1:])]
    manh3_plen[1][0] += 1

    # Dynamic programming starts
    for i in range(1, (m + 1)):
        for j in range(1, (n + 1)):
            tmp_seq = [log_trans_dic[S1[i - 1][0]][S2[j - 1][0]] + log_link_p[1], log_equil_dict[S2[j - 1][0]] + log_link_p[2]]
            tmp = [tmp_seq[0] + log_S_epi["".join(S1[i - 1][1:])] + Transition_g_sum(S1[i - 1], S2[j - 1], log_trans_dic, log_S_epi, params['hypN']),
                   tmp_seq[1] + log_S1_epi["".join(S2[j - 1][1:])] + log_S2_epi["".join(S1[i - 1][1:])]]
            max_t = max(tmp)

            tmp0 = log_lamb_mu + log_equil_dict[S1[i - 1][0]] + log_link_p[0]
            tmp1 = log_lamb_mu + log_equil_dict[S1[i - 1][0]]
            tmp2 = log_equil_dict[S2[j - 1][0]] + log_lamb_beta

            ent0 = tmp0 + log_S2_epi["".join(S1[i - 1][1:])] + manh3[0][j]
            ent1 = tmp1 + max_t + manh3[0][j - 1]
            ent2 = tmp2 + log_S1_epi["".join(S2[j - 1][1:])] + manh2[1][j - 1]

            plen_v0 = manh3_plen[0][j] + 1
            plen_v1 = manh3_plen[0][j - 1] + 1
            plen_v2 = manh2_plen[1][j - 1] + 1

            max_v3, max_v2, index_v3, index_v2 = Maxima([ent0, ent1, ent2], [plen_v0, plen_v1, plen_v2])

            manh3[1][j] = max_v3
            manh2[1][j] = max_v2
            plen_ulist = [plen_v0, plen_v1, plen_v2]
            manh3_plen[1][j] = plen_ulist[index_v3]
            manh2_plen[1][j] = plen_ulist[index_v2]

            if params['align_path']:
                max_it = tmp.index(max_t)
                Backtrack_mat(index_v3, max_it, backtrack0, i, j, 0)

        manh3[0] = manh3[1]
        manh3[1] = [Na] * (n + 1)

        manh2[0] = manh2[1]
        manh2[1] = [Na] * (n + 1)

        manh3_plen[0] = manh3_plen[1]
        manh3_plen[1] = [0] * (n + 1)

        manh2_plen[0] = manh2_plen[1]
        manh2_plen[1] = [0] * (n + 1)
        if i < m:
            manh3[1][0] = manh3[0][0] + log_lamb_mu + log_equil_dict[S1[i][0]] + log_S2_epi["".join(S1[i][1:])] + log_link_p[0]
            manh3_plen[1][0] = manh3_plen[0][0] + 1

    norm_max_list = [m1 / m2 for m1, m2 in zip(manh3[0], manh3_plen[0])]

    S.L = max(norm_max_list)
    j_start = norm_max_list.index(S.L)
    S.loc = j_start
    S.averagedL = sum(norm_max_list) / float(len(norm_max_list))
    if params['all_prob']:
        S.prob = norm_max_list
    if params['align_path']:
        Recons_path(S, backtrack0, j_start)
    return S

def manhattanWrapper(arg):
    return Manhattan(arg['S'], arg['X'], arg['log_trans_dic'], arg['log_equil_dict'], arg['params'])

def prepareManhattanParams(x, S, log_trans_dic, log_equil_dict, params):
    jobs = []
    for entry in S:
        job = dict()
        job['X'] = copy.deepcopy(x)
        job['S'] = entry
        job['log_trans_dic'] = copy.deepcopy(log_trans_dic)
        job['log_equil_dict'] = copy.deepcopy(log_equil_dict)
        job['params'] = copy.deepcopy(params)
        jobs.append(job)
    return jobs

def Manhattan_obj(x, S, e_dict_norm, weights, log_equil_dict, p_num, params):
    '''Distribute region pairs to difference processes for parallel computing.

    Args:
        x: the parameter vector.
        S: the list of region pairs. Each element is a HomoRegion object.
        e_dict_norm: the dictionary of equilibrium probabilities (not adjusted by weight).
        weights: the weights vector.
        log_equil_dict: the dictionary of log equilibrium probabilities (adjusted by weight).
        p_num: maximal number of processes to run.
        params: additional parameters that need to be passed to the worker.

    Returns:
        None. The function will update the S list directly.
    '''
    if min(x) < 0:
        print("minus p", file=sys.stderr)
        return float("inf")
    Log_transition_dic = Trans_matrix(len(S[0].S1[0]) - 1, x, e_dict_norm, weights)
    jobs = prepareManhattanParams(x, S, Log_transition_dic, log_equil_dict, params)

    with get_context("spawn").Pool(p_num) as p:
        mp_queue = p.map(manhattanWrapper, jobs)
        S[:] = mp_queue
        p.close()
        p.join()
    return


def Main():
    args = ParseArg()
    S = ReadInput(args.Input)

    params = dict()
    params['hypN'] = 1

    p_num = args.process_num

    x, weights, equil_dict, log_equil_dict = ReadParameters(args.equil_file)
    equil_dict_norm = copy.deepcopy(equil_dict)
    Mod_equilibrium(equil_dict, log_equil_dict, weights)
    print(weights, file=sys.stderr)

    # global S_epi, log_S_epi
    params['S_epi'], params['log_S_epi'] = Epi_equilibrium(equil_dict, log_equil_dict, len(S[0].S1[0]) - 1, weights)

    params['all_prob'] = args.out_allvec
    params['align_path'] = args.align_path

    t0 = time()

    Manhattan_obj(x, S, equil_dict_norm, weights, log_equil_dict, p_num, params)
    t1 = time()
    print("Done mapping input sequences.", file=sys.stderr)
    if params['all_prob']:
        with open(args.output, "w") as fout, open(args.out_allvec, "w") as fout2:
            for s in S:
                print("\t".join([s.name, str(s.L), str(s.averagedL), str(s.loc)]), file=fout)
                print(",".join([str(f) for f in [s.name] + s.prob]), file=fout2)
    else:
        with open(args.output, "w") as fout:
            for s in S:
                print("\t".join([s.name, str(s.L), str(s.averagedL), str(s.loc)]), file=fout)

    if params['align_path']:
        with open(params['align_path'], "w") as fout2:
            for s in S:
                Print_path(s, fout2)

    print("Time:%f" % ((t1 - t0) / 60), file=sys.stderr)


if __name__ == "__main__":
    Main()
