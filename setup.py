from setuptools import setup

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    author='Robyn Thiessen-Bock',
    author_email='thiessenbock@nceas.ucsb.edu',
    name='pdgstaging',
    version='0.1.0',
    description='PDG Visualization staging pipeline',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/PermafrostDiscoveryGateway/viz-staging',
    packages=['pdgstaging'],
    install_requires=[
        'numpy >= 1.2, < 2.0',
        'pandas >= 1.4, < 2.0',
        'geopandas >= 0.10, < 1.0',
        'morecantile >= 3.1, < 4.0',
        'Rtree >= 0.9, < 1.0',
        'filelock >= 3.6, < 4.0'
    ],
    python_requires='>=3.9, <4.0',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ],
    license='Apache Software License 2.0',
)
