from setuptools import setup
from setuptools import find_packages

setup(
    name='dna2vec',
    version='1.0.0',
    description='dna2vec: Consistent vector representations of variable-length k-mers',
    author='Patrick Ng',
    author_email='pn.appdev@gmail.com',
    url='https://github.com/pnpnpn/dna2vec',
    license='MIT',
    install_requires=[
        'gensim>=0.13,<1.0',
        'logbook',
    ],
    packages=find_packages())
