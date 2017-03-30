import random
import string
import resource
import logbook
import arrow
import numpy as np
import os

def split_Xy(df, y_colname='label'):
    X = df.drop([y_colname], axis=1)
    y = df[y_colname]
    return (X, y)

def shuffle_tuple(tup, rng):
    lst = list(tup)
    rng.shuffle(lst)
    return tuple(lst)

def shuffle_dataframe(df, rng):
    """
    this does NOT do in-place shuffling
    """
    return df.reindex(rng.permutation(df.index))

def random_str(N):
    return ''.join(random.SystemRandom().choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for _ in range(N))

def memory_usage():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1E6

def estimate_bytes(filenames):
    return sum([os.stat(f).st_size for f in filenames])

def get_output_fileroot(dirpath, name, postfix):
    return '{}/{}-{}-{}-{}'.format(
        dirpath,
        name,
        arrow.utcnow().format('YYYYMMDD-HHmm'),
        postfix,
        random_str(3))
