# load: basic Docker image
FROM python

# Create a "work folder"
WORKDIR /app

# Copy project in this work folder
COPY . /app

# RUN: Install VibTools package itself
RUN pip install /app/

# RUN: Installation of testing dependencies
RUN python -m pip install -r /app/requirements/requirements_unittests.txt

# RUN: Installation of quality checker
RUN pip install pytest-cov \
    && pip install docstr-coverage \
    && pip install pylint
