import numpy as np
import os
from Bio.Seq import Seq

NUCLEOTIDES = 'ACGT'

class Tuple4:
    def __init__(self, pos1, pos2, neg1, neg2):
        self.pos1 = pos1
        self.pos2 = pos2
        self.neg1 = neg1
        self.neg2 = neg2

def determine_out_filename(output_dir, fileroot, mode, extension='txt'):
    return os.path.join(output_dir, '{}.{}.{}'.format(fileroot, mode, extension))

def create_tuple4(kmer1, kmer2, kmer1_neg, kmer2_neg):
    """
    all inputs are list of single nucleotides, e.g. ['A', 'A', 'C']
    """
    return Tuple4(
        ''.join(kmer1),
        ''.join(kmer2),
        ''.join(kmer1_neg),
        ''.join(kmer2_neg))

def insert_snippet(seq, snippet, idx):
    """
    idx: 0 <= idx <= len(seq)]
    """
    split1 = seq[:idx]
    split2 = seq[idx:]
    return split1 + snippet + split2

def pairwise_key(v1, v2):
    return '{}:{}'.format(v1, v2)

def rand_kmer(rng, k_low, k_high=None):
    """
    k_low and k_high are inclusive
    """
    if k_high is None:
        k_high = k_low
    k_len = rng.randint(k_low, k_high + 1)
    return ''.join([NUCLEOTIDES[x] for x in rng.randint(4, size=k_len)])

def rand_nt(rng):
    return NUCLEOTIDES[rng.randint(4)]

def generate_revcompl_pair(k_low, k_high=None, rng=None):
    # TODO make params k_high and rng be required
    if k_high is None:
        k_high = k_low
    if rng is None:
        rng = np.random
    kmer = rand_kmer(rng, k_low, k_high)
    return (kmer, revcompl(kmer))

def revcompl(kmer):
    return str(Seq(kmer).reverse_complement())

def generate_1nt_mutation_4tuple(rng, k_len):
    kmer1 = list(rand_kmer(rng, k_len, k_len))
    kmer2 = list(rand_kmer(rng, k_len, k_len))

    idx = rng.randint(len(kmer1))
    original_nt = kmer1[idx]
    mutate_nt = rand_nt(rng)

    kmer1_neg = list(kmer1)
    kmer1_neg[idx] = mutate_nt

    kmer2[idx] = mutate_nt
    kmer2_neg = list(kmer2)
    kmer2_neg[idx] = original_nt

    return create_tuple4(kmer1, kmer2, kmer1_neg, kmer2_neg)
