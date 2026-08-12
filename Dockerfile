FROM python:3.11-slim

# Install transitive dependencies and build tools
RUN apt-get update \
    && apt-get install -y git libspatialindex-dev libgdal-dev libproj-dev build-essential g++ \
    && rm -rf /var/lib/apt/lists/*

# Install uv
RUN pip install uv

# Install pdgstaging from GitHub repo using uv
RUN uv pip install --system git+https://github.com/PermafrostDiscoveryGateway/viz-staging.git@feature-tiling-k8s

WORKDIR /app

CMD ["python3"]