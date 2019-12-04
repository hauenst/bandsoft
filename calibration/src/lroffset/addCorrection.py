from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# Add iteration to original offset correction to get better quality
iter_tdc = {};
iter_ftdc = {};
with open("../../include/corr_lr_offsets.txt","rb") as g, open("../../include/lr_offsets.txt","rb") as f, open("../../include/improved_lr_offsets.txt","w") as h:
	for line in g:
		parse = line.strip().split("\t")
		ID = int(parse[1])*100 + int(parse[0])*10 + int(parse[2])
		iter_tdc[ID] = float(parse[3])
		iter_ftdc[ID] = float(parse[4])

	for line in f:
		parse = line.strip().split("\t")
		ID = int(parse[1])*100 + int(parse[0])*10 + int(parse[2])

		new_off_tdc = float(parse[3]) + iter_tdc[ID]
		new_off_ftdc = float(parse[4]) + iter_ftdc[ID]

		h.write( "%i\t%i\t%i\t%f\t%f\t0\t0\n" % (int(parse[0]) , int(parse[1]), int(parse[2]), new_off_tdc, new_off_ftdc ) )
