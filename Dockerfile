FROM ubuntu AS builder
WORKDIR /project
RUN apt-get update && \
    apt-get install -y bash cmake git vim gcc wget python3 xterm openssh-server nvidia-visual-profiler nvidia-nsight imagemagick \
    libmagickwand-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev
RUN wget -q https://github.com/GPUOpen-Tools/CodeXL/releases/download/v2.6/codexl_2.6-302_amd64.deb
RUN dpkg -i codexl_2.6-302_amd64.deb
RUN useradd -m chapter13
RUN echo "chapter13\n chapter13\n" > passwd chapter13
RUN echo export PATH=/opt/CodeXL_2.6-302:${PATH} >> /home/chapter13/.bash_profile

RUN git clone --recursive https://github.com/essentialsofparallelcomputing/Chapter13.git

RUN bash

ENTRYPOINT ["bash"]
