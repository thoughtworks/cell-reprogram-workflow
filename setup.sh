# clone tr to some folder
# cd, run 


mkdir ./deps
git clone --depth=1 https://gitlab.lcsb.uni.lu/CBG/transsynw.git ./deps/transsynw

wget https://github.com/Lan-lab/SIGNET/raw/main/SIGNET_Tutorial/SIGNET.py -O ./deps/SIGNET.py

# # curl -s 'https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv' >> trrust_rawdata_human.tsv
# echo "Setup completed"