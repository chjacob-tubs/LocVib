# load: basic Docker image
FROM python:3.11.4

# Create a "work folder"
WORKDIR /app

# Copy project in this work folder
COPY . /app

# RUN: Installation of dependencies and package itself  
RUN  python -m pip install -r /app/requirements/requirements_unittests.txt \
     && pip install pytest-cov \
     && pip install /app/
