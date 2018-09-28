#!/usr/bin/python
import sys
#should be sequence centered at the peak summit
#check the middle part (100-200) if there are significant presence of repeats (N or a,g,c,t)
#if so then flag this sequence to be bad sequence

f = open(sys.argv[1])
tt = []
seq = {}
for l in f:
    l = l.rstrip("\n")
    if l.startswith(">"):
        tit = l[1:]
        seq[tit] = ""
        tt.append(tit)
        continue
    seq[tit] += l
f.close()

for tit in seq:
	seq[tit].replace(" ", "")

seq_len = int(sys.argv[2]) #sequence length
bad_s = set([])

st = int(seq_len / 3.0)
end = int(seq_len / 3.0 * 2.0)

for s in seq:
	num_N=0
	for i in range(st, end):
		#if seq[s][i]=="N":
		if seq[s][i]=="N" or \
		seq[s][i]=="a" or seq[s][i]=="g" or seq[s][i]=="c" or seq[s][i]=="t":
			num_N+=1
	if float(num_N) > 30: #originally 30
		bad_s.add(s)

for s in tt:
	if s in bad_s:
		tx = s.split(":")
		tx2 = tx[1].split("-")
		sys.stdout.write("%s\t%s\t%s\n" % (tx[0], tx2[0], tx2[1]))
