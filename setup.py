from setuptools import setup

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='pdgstaging',
    version='0.1.0',
    description='PDG Visualization staging pipeline',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/PermafrostDiscoveryGateway/viz-staging',
    py_modules=['pdgstaging'],
    install_requires=[
        'numpy==1.22.2',
        'pandas==1.4.1',
        'geopandas==0.10.2',
        'morecantile==3.1.0',
        'Rtree==0.9.7'
    ],
    python_requires='>=3.9',
)
