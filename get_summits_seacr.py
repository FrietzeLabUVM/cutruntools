#!/usr/bin/python

import sys
import os
import re

f = open(sys.argv[1])
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	coord = ll[-1]
	c1 = coord.split(":")[0]
	cx = coord.split(":")[1]
	t1 = int(cx.split("-")[0])
	if cx.split("-")[1]=="":
		t2 = t1
	else:
		t2 = int(cx.split("-")[1])
	
	mid = (t1+t2)/2
	sys.stdout.write("%s\t%d\t%d\n" % (c1, mid, mid+1))
f.close()
