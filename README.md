# ISOETES1 - A pipeline for classification of short-read plant metagenomic data from environmental or ancient DNA

## ISOETES1 is a copyright of Robert S. Harbert, 2017 and is shared according to the CC-BY-NC-SA license.

Code in this repository is associated with analyses done for:
Harbert, R.S. 2018. Algorithms and strategies in short-read shotgun 
	metagenomic reconstruction of plant communities. Applications in Plant Sciences [In Press].

NOTE: ANY FUTURE DEVELOPMENT OF THIS CODE WILL GO INTO A NEW REPOSITORY: ISOETES2. 

#This is a legacy pipeline for publication documentation ONLY, and is not meant to be suitable for reuse "as-is".
	
In read classification scripts QC and sga preprocessing steps were adapted from the HOLI pipeline (https://github.com/ancient-eDNA/Holi)

# Install dependencies:
Note: You will also need to make sure you have MegaBLAST installed.

```
mkdir ~/bin #local user bin. Could install to default locations if you have root privileges.

#in the ISOETES1 directory:
git clone https://github.com/ancient-eDNA/Holi
git clone https://github.com/DerrickWood/kraken
git clone https://github.com/HadrienG/Kraken_db_install_scripts.git #for updated database scripts
git clone https://github.com/infphilo/centrifuge.git
git clone https://github.com/marbl/Krona.git

#Get Jellyfish 1.*
wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz
tar -xzf jellyfish-1.1.11.tar.gz
rm jellyfish-1.1.11.tar.gz

#To install jellyfish in a custom directory:
cd jellfyish-1.1.11
./configure --prefix=$HOME/bin/jellyfish
make
make install

export PATH=$HOME/bin/jellyfish/bin:$PATH
export LD_LIBRARY_PATH=$HOME/bin/jellyfish/lib:$LD_LIBRARY_PATH
export MANPATH=$HOME/bin/jellyfish/share/man:$MANPATH
export PKG_CONFIG_PATH=$HOME/bin/jellyfish/lib/pkgconfig:$PKG_CONFIG_PATH
##OR APPEND ~/.bashrc

mv ./Kraken_db_install_scripts/kraken_shell_scripts/* ./kraken/scripts

##Kraken Installation:
mkdir kraken_bin
cd kraken
./install_kraken.sh ./kraken_bin

##Move kraken scripts to bin...
cp ./kraken_bin/* $HOME/bin/

cd Krona/KronaTools

./install.pl --prefix ./ --taxonomy ./taxonomy
./updateAccesions.sh

cd centrifuge
    make
    make install prefix=$HOME/bin



```



# Reference database
Get sequence data from NCBI and attach accession and taxon ID to each record in FASTA format
```
mkdir plastids_refseq 
#Note, the script below creates parallel download processes. This is not advised on wireless internet.
./download_referenc/ncbi2kraken.py -e 'your@email.address' -s 'chloroplast[All Fields] AND (refseq[filter] AND chloroplast[filter]' -o ./plastids_refseq
mkdir plants_genbank
./download_reference/ncbi2kraken.py -e 'your@email.address' -s 'genbank[filter] AND (plants[filter] AND biomol_genomic[PROP] AND ("200"[SLEN] : "5000"[SLEN]))' -o ./plants_genbank

#Merge files and convert to centrifuge compatible files
cat ./plastids_refseq/final.fa ./plants_genbank/final.fa > ref.fa

## Make centrifuge compatible file

sed "s/|kraken:taxid|/ /" ref.fa > centref.fa


#get accession to taxid map
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz
cut -f2,3 nucl_gb.accession2taxid > acc2tid.map


#Build Kraken kmer database (on 64 threads, adjust as needed)
kraken-build --add-to-library ./plastids_refseq/final.fa --db ./kdb_test
kraken-build --add-to-library ./plants_genbank/final.fa --db ./kdb_test
kraken-build --build --threads 64 --max-db-size 150 --db ./meta/kdb_test ##May require full path to database
kraken-build --clean --db ./kdb_test


#Build Centrifuge BWT index.
#using centref.fa
mkdir centrifuge ##May need to change name if you installed centrifuge in your ISOETES directory
##BEWARE centrifuge-build will require an enormous amount of RAM for building the database. Adjust bmax down if this is a problem.
centrifuge-download -o taxonomy taxonomy
centrifuge-build -p 16 --bmax 1619356040 --conversion-table acc2tid.map --taxonomy-tree ./taxonomy/nodes.dmp --name-table ./taxonomy/names.dmp ./centref.fa ./centrifuge/plants2



#Make blast databse:
cd blastdb
makeblastdb -in ./ref.fa -parse_seqids -dbtype nucl
cd..

```
# Metagenomic Simulations
```
#mkdir meta_sim
#subdirectories in meta_sim will be generated by the script at ./reproduce/meta_sim
#Short-read metagenomic simulated data will be generated by sampling the reference database and running wgsim as:
grep ">" ./centref.fa | awk -F " " '{print $2}' | sort | uniq | shuf | head -$nsample > sample.tids
cat sample.tids


sed 's/^/\\s/' sample.tids > sam.tids
sed 's/$/\\s/' sam.tids > samout.tids
awk -v RS='>' 'NR>1 {gsub("\n", ";", $0); sub(";$", "", $0); print ">"$0}' ./centref2.fa | grep -f samout.tids | tr ';' '\n' > ./sample.fa
rm sam.tids
rm samout.tids

cat sample.tids | sort
grep ">" ./sample.fa | awk -F " " '{print $2}' | sort | uniq

wgsim -N 100000 -e ./sample.fa ./out1.fq /dev/null > sim.log 
```


## Pedersen et al. (2016) data
```
##Get soil metagenomic sequences
mkdir fastq
cd fastq
wget $(cat ../urlist)
gunzip -f *.gz
cd .. 

```
## Metagenomic classification

### Simulation
```
mkdir meta_sim
#Run script reproduce/meta_sim as:
#will sample 10 taxa, classify, and output result to ./meta_sim/10_1
./reproduce/meta_sim simfile meta_sim/10_1 /full/path/to/ISOETES1 10

```

### Soil metagenomes from the Ice-Free Corridor
```
#For each file in fastq run:
./reproduce/kraken_search_fastq ./fastq/sample_ID.fastq soils.res /full/path/to/ISOETES1

```

## Look at compile.R and process_taxa.R in the reproduce directory for code to generate summary figures in R
