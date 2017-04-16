"""Setuptools entry point."""
import codecs
import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: MIT License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Software Development :: Libraries :: Python Modules'
]

dirname = os.path.dirname(__file__)

long_description = (
    codecs.open(os.path.join(dirname, 'README.md'), encoding='utf-8').read()
)

setup(
    name='dna2vec',
    version='1.0.0',
    description='dna2vec',
    long_description=long_description,
    author='Patrick Ng',
    author_email='pn.appdev@gmail.com',
    url='https://github.com/pnpnpn/dna2vec',
    packages=['dna2vec'],
    install_requires=[
        'gensim>=0.13,<1.0'
    ],
    classifiers=CLASSIFIERS)
