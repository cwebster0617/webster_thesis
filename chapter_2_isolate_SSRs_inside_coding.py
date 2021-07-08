from Bio import Seq
from Bio.Seq import Seq
import re, sys, argparse, os

#Arguments that can be altered in the terminal
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta_file_name', type=str)
parser.add_argument('-s', '--SciRoKo_output_file_name', type=str) #td file type
parser.add_argument('-l', '--SSR_minimum_repeat', type=int) #2=dinucleotide
parser.add_argument('-out_file', '--output_file_name', type=str) #Output file contains transcript name, motif, and positions only where SSRs occur in the coding regions
parser.add_argument('-out_fasta', '--output_fasta_file_names', type=str) #fasta files contain full gene sequence with SSR sequences that occur in the coding regions

args = parser.parse_args()

#FASTA file input
if args.fasta_file_name:
    fastaIn = args.fasta_file_name
else:
    print('Error: Need file containing fasta sequences')
    
#SciRoKo output file
if args.SciRoKo_output_file_name:
    SciRoKoIn = args.SciRoKo_output_file_name
else:
    print('Error: Need file containing SSR locations')
    
#SSR minimum repeat type
if args.SSR_minimum_repeat:
    Minrep = args.SSR_minimum_repeat
else:
    SSR_minimum_repeat = 2
    
#Output file with transcript name, SSR, start, stop
if args.output_file_name:
    FileOut = args.output_file_name
else:
    print('Error: Output file name not given, or already exists')

#FASTA output files, to output directory
if args.output_fasta_file_names:
    fastaOut = args.output_fasta_file_names
    if os.path.exists(fastaOut):
        print('output folder already exists, choose new folder')
        quit()
    else:
        os.mkdir(fastaOut)
else:
    print('Error: Output directory name not given')

ssr_min = Minrep



#Begin script
#PART 1: reading in fasta, translating, and storing info

f = open(fastaIn,'r')
lines = f.readlines()
f.close()

#Dictionary with transcript and coding region information
stored = {}
#transcript:[coding, str(frame), str(start), str(stop)]

transcript = []
seqs = []
#to store transcript name (start, stop, increments)
for l in range(0,len(lines),2):
    transcript.append(lines[l].strip())
    #defining 'key', transcript-1 gives first element in transcript
    key = transcript[-1].strip('>')
    #to store sequences (transcript+1 = sequences)
    s = lines[l+1].strip()
    seqs.append(s)
    #begin translation code
    translated = []
    for x in range(0,len(lines)):
        coding_dna = Seq(s[x:])
        translated.append(coding_dna.translate())
    coding = ''
    start = 0
    frame = 0
    stop = 0
    for t in range(0,len(lines)):
        c = translated[t].split('*')
        a = 0
        for r in c:
            if len(r) > len(coding):
                #if r is greater than coding, update coding with new r
                coding = r
                #if r is greater than coding, store t(translated) in frame
                frame = t
                #if r is greater than coding, multiply amino acids by 3 and add one for the stop codon that was split on, need nucleotide position
                start = a*3 + t
                #if r is greater than coding, use start and the length of coding x3, need nucleotide position
                stop = start + len(coding)*3
                #for dictionary (stored)
            a += len(r)+1
    #for dictionary (stored)
    stored[key] = [coding, int(frame), int(start), int(stop)]



#PART 2: reading in SSR data (transcript, start, stop positions) to tell if it falls in coding region

#from SciRoKo td output file
#this is where script can be altered to read regular BED file
m = open(SciRoKoIn,'r')
lines2 = m.readlines()
m.close()

info = ''.join(lines2).split('\r\n')

ssr1 = []

ssr_in_coding = []
ssr_out_of_coding = []

for i in info[1:]:
    sl = i.split('\t')
    print('test',sl)
    if len(sl) > 5 and len(sl[1]) >= ssr_min:
        #storing information from SciRoKo td output file (transcript, motif, start, stop)
        ssr1.append([sl[0],sl[1],int(sl[3]),int(sl[4])])

for s in ssr1:
    if s[0] in stored.keys():
        if s[2] > stored[s[0]][2] and s[3] < stored[s[0]][3]:
            ssr_in_coding.append(s)
        else:
            ssr_out_of_coding.append(s)



#PART 3: creating file output

forfasta = {}

#Writing file with info on just SSRs in coding region
c = open(FileOut, "w+")

c.write('Transcript'+'\t'+'SSR'+'\t'+'Start'+'\t'+'Stop'+'\n')
for r in ssr_in_coding:
    c.write('\t'.join(map(str,r))+'\n')

    
    
#PART 4: FASTA file outputs (to output directory)

#These fasta files can be used with Clustal Omega to visually see where SSRs occur
#Writing individual fasta files for transcripts (cycle through the forfasta key)

    if r[0] not in forfasta.keys():
        print(r[0])
        forfasta[r[0]] = ['>'+r[0],seqs[transcript.index('>'+r[0])]]
    
    forfasta[r[0]].append('>'+r[0])
    forfasta[r[0]].append(seqs[transcript.index('>'+r[0])][r[2]-1:r[3]])
c.close()

for k in forfasta.keys():
    y = open("%s/RTS_perfect_SSR_repeat_inside_coding_%s.txt"%(fastaOut,k), "w+")
    y.write('\n'.join(forfasta[k]))
    y.close
