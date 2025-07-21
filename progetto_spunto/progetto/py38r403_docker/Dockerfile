
ARG UBUNTU_RELEASE=20.04
FROM ubuntu:$UBUNTU_RELEASE
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ENV CRAN_MIRROR=https://cloud.r-project.org \
    CRAN_MIRROR_TAG=-cran40

USER root

COPY install_apt.sh    /opt/
# COPY install_pip.sh    /opt/
COPY install_rpacks.sh /opt/

RUN sh /opt/install_apt.sh 

RUN sh /opt/install_rpacks.sh       

RUN python -m venv /opt/python3_env 
RUN source /opt/python3_env/bin/activate 

### PIP INSTALL ###
RUN python -m pip --no-cache-dir install pip        --upgrade
RUN python -m pip --no-cache-dir install setuptools --upgrade
RUN python -m pip --no-cache-dir install wheel      --upgrade

RUN pip install --no-cache-dir jinja2
RUN pip install --no-cache-dir pytest
RUN pip install --no-cache-dir sphinx
RUN pip install --no-cache-dir tzlocal

RUN pip install --no-cache-dir torch>=1.8.0
# RUN pip install --no-cache-dir cudnn>=10.2
RUN pip install --no-cache-dir numpy==1.22.3
RUN pip install --no-cache-dir scanpy==1.9.1
RUN pip install --no-cache-dir squidpy
RUN pip install --no-cache-dir anndata==0.8.0
RUN pip install --no-cache-dir Pyrebase4
RUN pip install --no-cache-dir rpy2==3.4.1
RUN pip install --no-cache-dir pandas==1.4.2
RUN pip install --no-cache-dir scipy==1.8.1
RUN pip install --no-cache-dir scikit-learn==1.1.1
RUN pip install --no-cache-dir tqdm==4.64.0
RUN pip install --no-cache-dir matplotlib==3.4.2
RUN pip install --no-cache-dir jupyter
RUN pip install --no-cache-dir notebook

RUN pip install --no-cache-dir POT
RUN pip install --no-cache-dir GraphST

RUN pip install --no-cache-dir scikit-misc

RUN rm -rf /root/.cache
### PIP INSTALL ###

RUN rm -rf /tmp/* 
RUN apt-get remove --purge -y $BUILDDEPS 
RUN apt-get autoremove -y 
RUN apt-get autoclean -y 
RUN rm -rf /var/lib/apt/lists/*

EXPOSE 8888
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]

# docker run -it -p 8888:8888 -v .\shared:/shared py38r403