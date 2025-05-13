## Primer design for site saturation mutagenesis of cEF hand of STIM1
import Bio
from Bio.Seq import Seq
from datetime import datetime

## Sequence before cEF hand.
before = "gcagattataccgcaactacacgccacc".upper()


## Sequence of ORF to mutagenize
stim_cEF = "atggcccggcgctgcggccccgtggcgctgctccttggcttcggcctcctccggctgtgctcaggggtgtggggtggagaacagaaattgatctccgaagaagacctgtacgacgatgac".upper()

## Sequence after cEF hand.
after = "gataagggtacggatacagaggagcgg".upper()

## Degenerate codon of choice
deg = "NNK"

## Concatenate all of the above
fullstring = Seq(before + stim_cEF + after)
#print(fullstring)

## Numbers for For Loop
startnum = 0

## Number of residues * 3 + 3 (since not inclusive)
endnum = 40 * 3 + 3

plist = list()
plistf = list()
plistrc = list()
text_file = open("SSM_primer_EFhand_Output.tsv", "w")
before_offset = len(before)

for x in range(startnum, endnum, 3):
	pos = fullstring[(before_offset + x):(before_offset + x + 3)]
	plist.append(str(pos))
	# posf adds two strings - one before NNK and one after the 3 letters NNK replaces
	posf = fullstring[(before_offset + x - 18):(before_offset + x)] + deg + fullstring[(before_offset + x + 3):(before_offset + x + 3 + 12)]
	plistf.append(str(posf))
	posr = fullstring[(before_offset + x - 24):(before_offset + x)]
	posrc = posr.reverse_complement()
	plistrc.append(str(posrc))

ntcount = 0

##forward primers
for x in range(0, 40):
	print("CHRNB2_f_%s\t%s" % (x+1, plistf[x]))
	text_file.write("CHRNB2_f_%s\t%s\n" % (x+1, plistf[x]))
	ntcount = ntcount + len(plistf[x])

##reverse primers
for x in range(0, 40):
	print("CHRNB2_r_%s\t%s" % (x+2, plistrc[x]))
	text_file.write("CHRNB2_r_%s\t%s\n" % (x+1, plistrc[x]))
	ntcount = ntcount + len(plistrc[x])

print("total cost: $ %d" % (ntcount*0.1968))
text_file.close()