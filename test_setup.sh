

# dependencies for TransSynW
python -m pip install --upgrade pip
pip install pytest
pip install -r requirements.txt
# R -e 'install.packages(c("purrr","Rcpp","reshape2","gtools","tibble","stringr","stringili"), repos = "http://cran.us.r-project.org")'

# R -e 'install.packages(https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-arm64/contrib/4.3/RcisTarget_1.22.0.tgz, repos = NULL, type="source")'

# R -e 'install.packages(https://cran.r-project.org/src/contrib/tibble_3.2.1.tar.gz, repos = NULL, type="source")'
R 
install.packages('RcisTarget', repos='http://cran.us.r-project.org', dependencies=TRUE)
install.packages('dplyr', repos='http://cran.us.r-project.org', dependencies=TRUE)
install.packages('tibble', repos='http://cran.us.r-project.org', dependencies=TRUE)
install.packages('Rcpp','gtools','stingr','purrr', repos='http://cran.us.r-project.org', dependencies=TRUE)
# sudo apt install libboost-all-dev

# export CPATH=$CPATH:'/opt/homebrew/include'
# export LIBRARY_PATH=$LIBRARY_PATH:'/opt/homebrew/lib/'

cd ./data
curl -O https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/ 
cd ..

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