python -m pip --no-cache-dir install pip        --upgrade
python -m pip --no-cache-dir install setuptools --upgrade
python -m pip --no-cache-dir install wheel      --upgrade

pip install --no-cache-dir jinja2
pip install --no-cache-dir pytest
pip install --no-cache-dir sphinx
pip install --no-cache-dir tzlocal

pip install --no-cache-dir torch>=1.8.0
pip install --no-cache-dir cudnn>=10.2
pip install --no-cache-dir numpy==1.22.3
pip install --no-cache-dir scanpy==1.9.1
pip install --no-cache-dir squidpy
pip install --no-cache-dir anndata==0.8.0
pip install --no-cache-dir rpy2==3.4.1
pip install --no-cache-dir pandas==1.4.2
pip install --no-cache-dir scipy==1.8.1
pip install --no-cache-dir scikit-learn==1.1.1
pip install --no-cache-dir tqdm==4.64.0
pip install --no-cache-dir matplotlib==3.4.2
pip install --no-cache-dir jupyter
pip install --no-cache-dir notebook

rm -rf /root/.cache