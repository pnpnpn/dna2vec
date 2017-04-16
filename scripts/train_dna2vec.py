#!/usr/bin/env python3

import sys
sys.path.extend(['.', '..'])

import glob
import logbook
from logbook.compat import redirect_logging
import configargparse
import numpy as np
from Bio import SeqIO
from attic_util.time_benchmark import Benchmark
from attic_util import util
from attic_util.tee import Tee
from dna2vec.histogram import Histogram
from dna2vec.generators import SeqGenerator, KmerSeqIterable, SeqMapper, SeqFragmenter
from dna2vec.generators import DisjointKmerFragmenter, SlidingKmerFragmenter

from gensim.models import word2vec

class InvalidArgException(Exception):
    pass

class Learner:
    def __init__(self, out_fileroot, context_halfsize, gensim_iters, vec_dim):
        self.logger = logbook.Logger(self.__class__.__name__)
        assert(word2vec.FAST_VERSION >= 0)
        self.logger.info('word2vec.FAST_VERSION (should be >= 0): {}'.format(word2vec.FAST_VERSION))
        self.model = None
        self.out_fileroot = out_fileroot
        self.context_halfsize = context_halfsize
        self.gensim_iters = gensim_iters
        self.use_skipgram = 1
        self.vec_dim = vec_dim

        self.logger.info('Context window half size: {}'.format(self.context_halfsize))
        self.logger.info('Use skipgram: {}'.format(self.use_skipgram))
        self.logger.info('gensim_iters: {}'.format(self.gensim_iters))
        self.logger.info('vec_dim: {}'.format(self.vec_dim))

    def train(self, kmer_seq_generator):
        self.model = word2vec.Word2Vec(
            sentences=kmer_seq_generator,
            size=self.vec_dim,
            window=self.context_halfsize,
            min_count=5,
            workers=4,
            sg=self.use_skipgram,
            iter=self.gensim_iters)

        # self.logger.info(model.vocab)

    def write_vec(self):
        out_filename = '{}.w2v'.format(self.out_fileroot)
        self.model.save_word2vec_format(out_filename, binary=False)

def run_main(args, inputs, out_fileroot):
    logbook.info(' '.join(sys.argv))
    if not args.debug:
        import logging
        logging.getLogger('gensim.models.word2vec').setLevel(logging.INFO)

    np.random.seed(args.rseed)

    benchmark = Benchmark()

    if args.kmer_fragmenter == 'disjoint':
        kmer_fragmenter = DisjointKmerFragmenter(args.k_low, args.k_high)
    elif args.kmer_fragmenter == 'sliding':
        kmer_fragmenter = SlidingKmerFragmenter(args.k_low, args.k_high)
    else:
        raise InvalidArgException('Invalid kmer fragmenter: {}'.format(args.kmer_fragmenter))

    logbook.info('kmer fragmenter: {}'.format(args.kmer_fragmenter))

    histogram = Histogram()
    kmer_seq_iterable = KmerSeqIterable(
        args.rseed_trainset,
        SeqGenerator(inputs, args.epochs),
        SeqMapper(),
        SeqFragmenter(),
        kmer_fragmenter,
        histogram,
    )

    learner = Learner(out_fileroot, args.context, args.gensim_iters, args.vec_dim)
    learner.train(kmer_seq_iterable)
    learner.write_vec()

    histogram.print_stat(sys.stdout)

    benchmark.print_time()

def main():
    argp = configargparse.get_argument_parser()
    argp.add('-c', is_config_file=True, help='config file path')
    argp.add_argument('--kmer-fragmenter', help='disjoint or sliding', choices=['disjoint', 'sliding'], default='sliding')
    argp.add_argument('--vec-dim', help='vector dimension', type=int, default=12)
    argp.add_argument('--rseed', help='general np.random seed', type=int, default=7)
    argp.add_argument('--rseed-trainset', help='random seed for generating training data', type=int, default=123)
    argp.add_argument('--inputs', help='FASTA files', nargs='+', required=True)
    argp.add_argument('--k-low', help='k-mer start range (inclusive)', type=int, default=5)
    argp.add_argument('--k-high', help='k-mer end range (inclusive)', type=int, default=5)
    argp.add_argument('--context', help='half size of context window (the total size is 2*c+1)', type=int, default=4)
    argp.add_argument('--epochs', help='number of epochs', type=int, default=1)
    argp.add_argument('--gensim-iters', help="gensim's internal iterations", type=int, default=1)
    argp.add_argument('--out-dir', help="output directory", default='../dataset/dna2vec/results')
    argp.add_argument('--debug', help='', action='store_true')
    args = argp.parse_args()

    if args.debug:
        out_dir = '/tmp'
        log_level = 'DEBUG'
    else:
        out_dir = args.out_dir
        log_level = 'INFO'

    inputs = []
    for s in args.inputs:
        inputs.extend(list(glob.glob(s)))

    mbytes = util.estimate_bytes(inputs) // (10 ** 6)
    out_fileroot = util.get_output_fileroot(
        out_dir,
        'dna2vec',
        'k{}to{}-{}d-{}c-{}Mbp-{}'.format(
            args.k_low,
            args.k_high,
            args.vec_dim,
            args.context,
            mbytes * args.epochs,  # total Mb including epochs
            args.kmer_fragmenter))

    out_txt_filename = '{}.txt'.format(out_fileroot)
    with open(out_txt_filename, 'w') as summary_fptr:
        with Tee(summary_fptr):
            logbook.StreamHandler(sys.stdout, level=log_level).push_application()
            redirect_logging()
            run_main(args, inputs, out_fileroot)

if __name__ == '__main__':
    main()
