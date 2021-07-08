import re,random, os, sys, argparse

#Arguments that can be altered in the terminal
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta_file_name', type=str) #look up add.argument = required
parser.add_argument('-v', '--vcf_file_name', type=str) #look up add.argument = required. Need to have in order to run.
parser.add_argument('-ad', '--allele_depth', type=int)
parser.add_argument('-sd', '--sample_depth', type=int)
parser.add_argument('-o', '--output_directory', type=str) #send files to an optput directory

args = parser.parse_args()

#fasta file
if args.fasta_file_name:
	fastaIn = args.fasta_file_name
else:
	print('Error: Need file containing fasta sequences')

#vcf file
if args.vcf_file_name:
	vcfIn = args.vcf_file_name
else:
	print('Error: Need file containing variant call format')

#allele depth
if args.allele_depth:
	allele_depth = args.allele_depth
else:
	allele_depth = 60

#sample depth
if args.sample_depth:
	sample_depth = args.sample_depth
else:
	sample_depth = 8

#output file make directory
if args.output_directory:
    Outdir = args.output_directory
    if os.path.exists(Outdir):
        print('output folder already exists, choose new folder')
        quit()
    else:
        os.mkdir(Outdir)
else:
	print('Error: Output directory name not given')



#PART 1
#Bring in sequence dictionary (find_seq.py)

f = open(fastaIn, 'r')
lines = f.readlines()
f.close()

#cleaning up sequences
seq = ''.join(lines)[4:].split('>')

search = {}

#storing headers in dictionary
for s in seq:
	temp = s.split('\n')
	head = temp[0]
	seqs = temp[1]
	search[head] = seqs



#PART 2
#Reading in VCF and determining if allele and sample depth qualify per parameters
positions = []

g = ''
#Reading in VCF
#For loop to go through sites on protein
with open(vcfIn, 'r') as fp:
    l = fp.readline()
    while l[0:2] == '##':
        l = fp.readline()
    #headers of individual sequence locations
    header2 = l.split('\t')[9:]

    head =[]
    all = []
    #creating headers for transformed sequence data
    #'a' for reference and 'b' for alternate mutations
    for h in header2:
        h = h.split('.')[0]
        head.append('>' + h + 'a')
        head.append('>' + h + 'b')

#Begin code for AD and SD parameters
    l = fp.readline()
    while l:
        sl = l.strip().split('\t')
        
        #determine AD
        if sl[0] in search.keys():
            AD1 = sl[7].split(';')
            AD2 = AD1[1].split('=')
            ADfinal = AD2[1].split(',')
            ADtotal = 0
            for x in range(0,len(ADfinal)):
                ADtotal += int(ADfinal[x])
            if ADtotal >= allele_depth and sl[0] in search.keys():
            
                #if the allele depth qualifies, then move to sample depth
                #determine SD
                SD1 = sl[9:]
                SD = 0
                for s in SD1:
                    if re.search('[1-9]',s):
                        SD += 1
                if SD >= sample_depth:
                #if both allele depth and sample depth qualify, transcript name will be printed in the terminal
                    print(sl[0])
                    if sl[0] != g:
                        #Wrap up old
                        if g != '':
                        
                            #writing individual output files with transformed data
                            f = open("%s/TRINITY_%s_individ_sequences.txt"%(Outdir,g),"w+")
                            for a in range (0,16):
                                f.write(head[a] + '\n')
                                f.write(''.join(all[a]) + '\n')
                            f.close()
                        #Start new
                        g=sl[0]
                        all = []
                        print(search[g])
                        for r in range(0,16):
                            all.append(list(search[g]))
                    print(all.append(list(search[g])))



#PART 3
#Begin code for data transformation, make 16 copies of each sequence, 2 for each sample (reference and alternate)
                    #Get info for all samples
                    ref = l.split('\t')[3]
                    alt = l.split('\t')[4].split(',')[0]

                    #For loop to go through samples
                    i=0
                    pos = int(l.split('\t')[1])-1
                    
                    #Begin boolean type for SNP BED file
                    Boo = False
                    for t in l.split('\t')[9:]:
                        st = t.split(':')[-1].split(',')
                        #assigning variables for reference and alternate nucleotides
                        trc = int(st[0])
                        tac = int(st[1])
                        
                        #for heterozygotes, randomly assigned to use ref or alt
                        if trc != 0 and tac != 0:
                            if float(trc) / (trc+tac) > 0.4 and float(trc) /(trc+tac) < 0.7: #hetero
                                Boo = True
                                r = random.random()
                                if r >= 0.5:
                                    all[i][pos] = alt #alt
                                else:
                                    all[i+1][pos] = alt #ref
                                    
                            #homozygous for the alternate
                            elif float(tac) / trc > 0.2: #homo alt
                                all[i][pos] = alt
                                all[i+1][pos] = alt
                                Boo = True
                        elif trc == 0 and tac != 0: #homo alt
                            all[i][pos] = alt
                            all[i+1][pos] = alt
                            Boo = True
                        i+=2
                        
                    #Ending boolean loop and storing in "positions" variable
                    if Boo == True:
                        positions.append(sl[0]+'\t'+l.split('\t')[1])
        l = fp.readline()



#PART 4
#Writing positions of SNPs to BED file
#BED file contains positions of SNPs in transcript and how many SNPs each transcript contains that meet the parameters
#Boolean type parameters, True = hetero and homo for the alternate, no homo for the reference
v = open("SNP_positions.txt", "w+")
v.write('\n'.join(positions))
v.close()



#PART 5
#Write to Output directory (with individual fasta files)
if g != '':
	f = open("%s/TRINITY_%s_individ_sequences.txt"%(Outdir,g),"w+")
	for a in range (0,16):
		f.write(head[a] + '\n')
		f.write(''.join(all[a]) + '\n')
	f.close()

