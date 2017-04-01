# dna2vec

**Dna2vec** is an open-source library to train distributed representations
of variable-length k-mers.

For more information, please refer to the paper: [dna2vec: Consistent vector representations of variable-length k-mers](https://arxiv.org/abs/1701.06279)

Installation
---

Note that this implementation has only been tested on Python 3.5.3, but we welcome any
contributions and bug reporting to make it more accessible.

1. Clone the `dna2vec` repository: `git clone https://github.com/pnpnpn/dna2vec`
2. Install Python dependencies: `pip3 install -r requirements.txt`
3. Test the installation: `python3 ./scripts/train_dna2vec.py -c configs/small_example.yml`

Training dna2vec embeddings
---

1. Download `hg38` from <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz>.
    This will take a while as it's 938MB.
2. Untar with `tar -zxvf hg38.chromFa.tar.gz`. You should see at least 22
    FASTA files `chr1.fa`, `chr2.fa`, ..., `chr22.fa`.


Reading pretrained dna2vec
---

You can read pretrained dna2vec vectors within `pretrained/dna2vec-*.w2v` using
the class `MultiKModel` in `dna2vec/multi_k_model.py`. For example:

```
from dna2vec.multi_k_model import MultiKModel

filepath = 'pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
mk_model = MultiKModel(filepath)
```

Contribute
---
I would love for you to fork and send me pull request for this project.
Please contribute.

License
---
This software is licensed under the [MIT license](http://en.wikipedia.org/wiki/MIT_License)
