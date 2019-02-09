#!/usr/bin/python
import sys
f = open(sys.argv[1])
outfile=sys.argv[2]
fw = open(outfile, "w")

freq = {}
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	c1,s1,e1,c2,s2,e2 = ll[:6]
	sign1,sign2 = ll[-2:]
	if c1==".":
		continue
	if sign1=="+":
		#frag ends are s1, e2
		fw.write("%s\t%d\t%d\n" % (c1, int(s1), int(s1)+1))
		fw.write("%s\t%d\t%d\n" % (c1, int(e2), int(e2)+1))

		freq.setdefault(int(s1), 0)
		freq[int(s1)]+=1
		freq.setdefault(int(e2), 0)
		freq[int(e2)]+=1
	else:
		#frag ends are s2, e1
		fw.write("%s\t%d\t%d\n" % (c1, int(s2), int(s2)+1))
		fw.write("%s\t%d\t%d\n" % (c1, int(e1), int(e1)+1))

		freq.setdefault(int(s2), 0)
		freq[int(s2)]+=1
		freq.setdefault(int(e1), 0)
		freq[int(e1)]+=1
f.close()
fw.close()

'''
for p in range(min(freq.keys()), max(freq.keys())+1):
	if p not in freq:
		#continue
		sys.stdout.write("%d\t%d\n" % (p, 0))
		continue
	sys.stdout.write("%d\t%d\n" % (p, freq[p]))
'''
