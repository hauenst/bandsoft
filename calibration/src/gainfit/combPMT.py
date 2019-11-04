from __future__ import division
import sys
import numpy as np


if len(sys.argv)!=3:
	print 'Invalid number of arguments. Please use:\n\tpython combPMT.py [OutputOfGainFit.txt] [FixedOutput.txt]'
	exit(-1)


Map = {}
with open(sys.argv[2],'w') as g:

	g.write("#Sector\tLayer\tComponent\tHV left\tHV right\n");
	with open(sys.argv[1],'rb') as a:

		for line in a:
			if '#' in line: continue

			li=line.strip().split(" ")
			
			hashCode = li[0]+li[1]+li[2]
			if( li[4] == 'inf' ): li[4] = 0.

			if hashCode in Map:
				Map[hashCode].append(int(float(li[4])))

			else:
				Map[hashCode] = [int(float(li[4]))]


	for key in Map:
                #print (key[0],key[1],key[2],Map[key][0])
                if int(key[1]) == 6:
			g.write("%s\t%s\t%s\t%d\n" % (key[0],key[1],key[2],Map[key][0]))
		else:
			g.write("%s\t%s\t%s\t%d\t%d\n" % (key[0],key[1],key[2],Map[key][0],Map[key][1]))
        
