#!/home/rharbert/miniconda2/bin/python
##A script to get data from NCBI, parse out taxon numbers and convert to kraken compatible fasta files.
## Writes all records as fasta to final.fa in working directory.


import sys, getopt

def main(argv):
	entrezemail = ''
	searchterm = ''
	outdir = ''
	try:
		opts, args = getopt.getopt(argv,"hs:e:o:",["search=","email=","outdir="])
	except getopt.GetoptError:
		print 'ncbi2kraken.py -s <ncbi_searchstring> -e <email> -o <outdir>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'ncbi2kraken.py -s <ncbi_searchstring> -e <email> -o <output_directory>'
			sys.exit()
		elif opt in ("-e", "--email"):
			entrezemail = arg
		elif opt in ("-s", "--search"):
			searchterm = arg
		elif opt in ("-o", "--odir"):
			outdir = arg

	return (entrezemail, searchterm, outdir)
		#print 'email is "', entrezemail
		#print 'searchterm is"', searchterm
		#print 'Output dir is "', outdir


if __name__ == "__main__":
   (email, search, outd) = main(sys.argv[1:])

import os
if not os.path.exists(outd):
	os.makedirs(outd)


from Bio import Entrez
Entrez.email=email
handle = Entrez.esearch(
	db="nucleotide", retmax=1000000, 
	term=search, 
	idtype="acc",
	usehistory="y")
record = Entrez.read(handle)
handle.close()



from joblib import Parallel, delayed
import multiprocessing

zz=200; #get in 200 record chunks
cuts = (len(record['IdList'])/zz)+1
inputs = range(0, cuts) 

def processInput(i):
	z = (i*zz)
	fetch = Entrez.efetch(db='nucleotide', id=",".join(record['IdList'][z:z+zz]), rettype='gb', retmode='text')
	fast = fetch.read()
	arfast = fast.split("//\n\n");
	seqs = [];
	taxid = [];
	for pi in arfast:
		if len(pi)>20:
			seqs.append(pi);
			arpi = pi.split("taxon:");
			nexpi = arpi[1].split("\"");
			taxid.append(nexpi[0]);
	return (seqs, taxid)

#num_cores = multiprocessing.cpu_count()
num_cores = 8;   
results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)
taxid = [];
seqs = [];
for n in range(0, len(results)):
	taxid = taxid+results[n][1]
	seqs = seqs+results[n][0]

fast = "//\n" .join(seqs)
targfile = outd + '/seqs.gb'
target = open(targfile, 'w')
target.write(fast)
target.close()

from Bio import SeqIO
in1 = outd + '/seqs.fasta'
count = SeqIO.convert(targfile, "genbank", in1, "fasta")
print("Converted %i records" % count)


infile = open(in1, 'r')
ffast = infile.read()
infile.close()
arffast = ffast.split("\n>")
n = 0;
arrout = [];
errout = [];
for line in arffast:
	if(taxid[n].isdigit()):
		##Don't add records that have not taxon ID attached
		str2 = "|kraken:taxid|%s" % taxid[n]
		str2 = str2 + " "
		ss= line.replace(" ", str2, 1)
		ss = ">" + ss
		arrout.append(ss); ##This array will be written to a fasta file compatible with Kraken DB building
		n=n+1
	else:
		str3 = "kraken:taxid|%s" % taxid[n]
		str3 = str3 + " "
		ss3 = line.replace(" ", str3, 1)
		ss3 = ">" + ss
		errout.append(ss3)
		n=n+1


outstr = "\n\n".join(arrout)
finale = outd + "/final.fa"
out = open(finale, 'w')
out.write(outstr)
out.close()	

outstr = "\n\n".join(errout)
finerr = outd + "/finerr.fa"
outerr = open(finerr, 'w')
outerr.write(outstr)
outerr.close()
exit()


