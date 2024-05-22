FROM python:3.9-alpine

# Install transitive dependencies
# RUN apk update && apk add git cmake make gcc
# RUN apk update && apk add git build-base
RUN apk update && apk add cmake alpine-sdk

# Install libspatialindex
# https://libspatialindex.org/en/latest/install.html
# TODO: Can we install with Conda? https://anaconda.org/conda-forge/libspatialindex
RUN git clone https://github.com/libspatialindex/libspatialindex.git
WORKDIR libspatialindex
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr
# RUN cmake # -DCMAKE_INSTALL_PREFIX=/usr
RUN make
RUN make install

# Install pdgstaging
RUN pip install git+https://github.com/PermafrostDiscoveryGateway/viz-staging.git
