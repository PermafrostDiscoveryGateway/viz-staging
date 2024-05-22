FROM python:3.9-alpine

# Install transitive dependencies
RUN apk update && apk add cmake alpine-sdk gdal proj

# Install libspatialindex
# https://libspatialindex.org/en/latest/install.html
# TODO: Can we install with Conda?
#       https://anaconda.org/conda-forge/libspatialindex
RUN git clone https://github.com/libspatialindex/libspatialindex.git
WORKDIR libspatialindex
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr \
 && make \
 && make install

# Install pdgstaging
# TODO: If we make pdg-staging conda-installable, we can install gdal, proj,
#       and libspatialindex that way too.
RUN pip install git+https://github.com/PermafrostDiscoveryGateway/viz-staging.git
