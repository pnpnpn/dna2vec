import pytest
from dna2vec.multi_k_model import MultiKModel

K_LOW = 3
K_HIGH = 8

kmer_samples = ['AAA', 'ACGT', 'ACGTACGT']
cosine_testdata = [
    ('AAA', 'AAAA', 'CCCC'),
    ('ACGA', 'ACGT', 'TTTT'),
    ('ACGT', 'ACGTAA', 'TTT'),
]

@pytest.fixture(scope="module")
def mk_model():
    filepath = 'pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
    return MultiKModel(filepath)

@pytest.mark.parametrize('k_len', list(range(K_LOW, K_HIGH + 1)))
def test_num_of_vectors(k_len, mk_model):
    assert len(mk_model.model(k_len).vocab) == 4 ** k_len

@pytest.mark.parametrize('kmer', kmer_samples)
def test_cosine_dist_to_itself(kmer, mk_model):
    assert abs(mk_model.cosine_distance(kmer, kmer) - 1.0) < 1e-10

@pytest.mark.parametrize('kmer0, kmer_greater, kmer_less', cosine_testdata)
def test_cosine_dist_sanity(kmer0, kmer_greater, kmer_less, mk_model):
    assert mk_model.cosine_distance(kmer0, kmer_greater) > mk_model.cosine_distance(kmer0, kmer_less)
