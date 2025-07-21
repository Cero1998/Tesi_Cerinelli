apt-get update --yes
apt-get upgrade --yes
apt install -y --no-install-recommends \
    software-properties-common         \
    dirmngr                            \
    lsb-release                        \
    wget

##### PYTHON 3.8 #####
apt-get update
apt install -y curl build-essential
add-apt-repository ppa:deadsnakes/ppa
apt install -y            \
  python3.8               \
  python3.8-dev           \
	python3.8-distutils     \
	python3.8-venv
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3.8 get-pip.py
ln -s /usr/bin/python3.8 /usr/bin/python
######################

##### R 4.0.3 #####
wget -qO- "${CRAN_MIRROR}"/bin/linux/ubuntu/marutter_pubkey.asc \
    | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

add-apt-repository \
        --yes \
        "deb ${CRAN_MIRROR}/bin/linux/ubuntu/ $(lsb_release -c -s)${CRAN_MIRROR_TAG}/"

apt-get update -qq
VERSION="4.0.3-1.2004.0"
apt install -y --no-install-recommends \
  r-base-core=${VERSION}               \
  r-base-html=${VERSION}               \
  r-doc-html=${VERSION}                \
  r-base-dev=${VERSION}
###################

rm -rf /var/lib/apt/lists/*