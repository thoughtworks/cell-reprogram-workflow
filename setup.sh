
# TransSynW dependencies

mkdir ./cmi/dependencies
cd ./cmi/dependencies
touch __init__.py
git clone --depth=1 https://gitlab.lcsb.uni.lu/CBG/transsynw.git transsynw

# Signet dependencies
wget https://github.com/Lan-lab/SIGNET/raw/main/SIGNET_Tutorial/SIGNET.py -O SIGNET.py

cd -
pip install .
