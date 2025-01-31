# Use the official Python 3.9 image from the Docker Hub
FROM python:3.11

# Set the working directory in the container
WORKDIR /sclive

# Copy the current directory contents into the container at /app
COPY . /sclive


RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="${PATH}:/root/.local/bin"
RUN ~/.local/bin/poetry install