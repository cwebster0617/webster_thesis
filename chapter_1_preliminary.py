import re,random, argparse

#looking for mutations in general across all Phytochrome genes

#Arguments that can be altered in the terminal
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--VCF_input_file', type=str, required=True) #input file
parser.add_argument('-o', '--parsed_output_file', type=str, required=True) #output file
parser.add_argument('-a', '--allele_depth', type=int, required=True)
parser.add_argument('-m', '--min_allele_freq', type=float, required=True)

args = parser.parse_args()

#VCF input file
if args.VCF_input_file:
    vcfIn = args.VCF_input_file
else:
    print('Error: Need VCF file')

#parsed output file (showing there are mutations in samples)
if args.parsed_output_file:
    fileOut = args.parsed_output_file
else:
    print('Error: Output file already exists')
    
#allele depth (only care about transcripts with _% for alternate)
if args.allele_depth:
    Depth = args.allele_depth
else:
    print('Error: No mutations at that allele depth')
    
#minimum allele frequency
if args.min_allele_freq:
    Freq = args.min_allele_freq
else:
    print('Error: No mutations at that allele frequency')

#Begin script
#PART 1 reading in and parsing FASTA file
f = open(vcfIn,'r')
lines = f.readlines()
f.close()
lines
temp = ''.join(lines).split('#')[-1]
geno = temp.strip().split('\n')

need = [] #list of sample lists
perc = [] #alternate numbers

#making a list of lists for samples
for g in geno[1:]:
	ge = g.strip().split('\t')
	need.append(ge)



#PART 2 creating output file with genes containing mutations
#Writing samples of interest (>= _allele depth and >_% for alt min allele freq)
f = open(fileOut,"w+")
f.write(geno[0] + '\n')
focus = []



#PART 3 calculating allele depth and minimum allele frequencey
#focus on 7th element where mutation quality is located
for ge in need:
	if ge[4] != '<*>':
		if ge[7].split(';')[0] == 'INDEL':
			temp = ge[7].split(';')[4].strip('AD=').split(',')
			if int(temp[0]) + int(temp[1]) >= Depth and float(temp[1]) / (int(temp[0]) + int(temp[1])) > Freq:
				focus.append([temp[0],temp[1]])
				f.write('\t'.join(ge) + '\n')
		else:
			temp = ge[7].strip().split(';')[1].strip('AD=').split(',')
			perc.append(temp)
			if int(temp[0]) + int(temp[1]) >= Depth and float(temp[1]) / (int(temp[0]) + int(temp[1])) > Freq:
				focus.append([temp[0],temp[1]])
				f.write('\t'.join(ge) + '\n')
f.close()
