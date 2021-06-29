# syntax=docker/dockerfile:1

FROM ubuntu:latest
WORKDIR /app
RUN apt-get update && apt-get install -y python3 python3-pip bowtie2
COPY requirements.txt requirements.txt
RUN ["pip3", "install", "-r", "requirements.txt"]
COPY extender extender
COPY build.sh build.sh
COPY setup.py setup.py
RUN ./build.sh
ENTRYPOINT ["python3", "extender/extender_wrapper.py"]
