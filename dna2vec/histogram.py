from collections import Counter
import logbook

class Histogram:
    def __init__(self):
        self.kmer_len_counter = Counter()
        self.nb_kmers = 0
        self.logger = logbook.Logger(self.__class__.__name__)

    def add(self, seq):
        """
        seq - array of k-mer string
        """
        for kmer in seq:
            self.kmer_len_counter[len(kmer)] += 1
            self.nb_kmers += 1

    def print_stat(self, fptr):
        for kmer_len in sorted(self.kmer_len_counter.keys()):
            self.logger.info('Percent of {:2d}-mers: {:3.1f}% ({})'.format(
                kmer_len,
                100.0 * self.kmer_len_counter[kmer_len] / self.nb_kmers,
                self.kmer_len_counter[kmer_len],
            ))

        total_bps = sum([l * c for l, c in self.kmer_len_counter.items()])
        self.logger.info('Number of base-pairs: {}'.format(total_bps))
