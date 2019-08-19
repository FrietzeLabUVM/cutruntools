#!/usr/bin/python

import sys

f = open(sys.argv[1])
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	ct = int(float(ll[-1]))
	if ct==0: continue
	sys.stdout.write("%s\t%s\t%s\t%d\n" % (ll[0], ll[1], ll[2], ct))
f.close()


