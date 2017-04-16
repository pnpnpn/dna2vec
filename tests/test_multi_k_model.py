import pytest
from dna2vec.multi_k_model import MultiKModel

K_LOW = 3
K_HIGH = 8

@pytest.fixture(scope="module")
def mk_model():
    filepath = 'pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
    return MultiKModel(filepath)

@pytest.mark.parametrize('k_len', list(range(K_LOW, K_HIGH + 1)))
def test_num_of_vectors(k_len, mk_model):
    assert len(mk_model.model(k_len).vocab) == 4 ** k_len
