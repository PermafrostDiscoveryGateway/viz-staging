FROM python:3.9-slim

# Install transitive dependencies
RUN apt-get update \
 && apt-get install -y git libspatialindex-dev libgdal-dev libproj-dev

# Install pdgstaging from GitHub repo
RUN pip install git+https://github.com/PermafrostDiscoveryGateway/viz-staging.git

WORKDIR /app

CMD ["python3"]
