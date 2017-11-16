############################################################
# Dockerfile to build DNAscan pipeline
############################################################


FROM ubuntu:latest
MAINTAINER Alfredo Iacoangeli "alfredo.iacoangeli@kcl.ac.uk"

RUN apt-get update
RUN apt-get install -y git
RUN cd /home/
RUN git clone https://github.com/snewhouse/DNA-NGS_scan.git
RUN cd DNA-NGS_scan
RUN bash scripts/install_dependencies_hg38_docker.sh /home/DNA-NGS_scan/local/ /home/DNA-NGS_scan/ 
RUN source ~/.bashrc

EXPOSE 8080

WORKDIR /home/DNA-NGS_scan/

ENTRYPOINT cat /home/DNA-NGS_scan/docker/welcome_message.txt
