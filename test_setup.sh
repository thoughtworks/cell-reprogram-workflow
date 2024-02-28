

# dependencies for TransSynW
python -m pip install --upgrade pip
pip install pytest
pip install -r requirements.txt
R -e 'install.packages(c("gtools","Matrix", "nibble","dplyr","stringr", "purrr","Rcpp","reshape2","umap", "pheatmap", "igraph","GGally","ggplot2", "RcisTarget","AUCell"), repos = "http://cran.us.r-project.org")'


mkdir ./craft/dependencies
cd ./craft/dependencies
touch __init__.py
git clone --depth=1 https://gitlab.lcsb.uni.lu/CBG/transsynw.git transsynw

# dependencies for Signet
wget https://github.com/Lan-lab/SIGNET/raw/main/SIGNET_Tutorial/SIGNET.py -O SIGNET.py

cd -

# if the user wants to download the latest data from TRRUST database
# # curl -s 'https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv' >> trrust_rawdata_human.tsv
# echo "Setup completed"