from distutils.core import setup, Extension, DEBUG


setup(name = 'vectorizer', version = '3.0',
    description = 'Python Package for Genome Vectorizer',
    ext_modules = [Extension('vectorizer', sources = ['vectorizer_bind.cpp','vectorizer.cpp'])]
    )