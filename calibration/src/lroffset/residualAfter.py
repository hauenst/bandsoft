from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# Add iteration to original offset correction to get better quality
iter_tdc = [];
iter_ftdc = [];
IDs = [];
with open("../../include/corr2nd_lr_offsets.txt","rb") as g:
	for line in g:
		parse = line.strip().split("\t")
		ID = int(parse[1])*100 + int(parse[0])*10 + int(parse[2])
		iter_tdc.append( float(parse[3]) )
		iter_ftdc.append( float(parse[4]) )
		IDs.append(ID)



plt.scatter(IDs,iter_tdc,color='red',label='TDC')
plt.scatter(IDs,iter_ftdc,color='blue',label='FADC')
plt.ylabel("Center of L-R Distribution [ns]",fontsize=16)
#plt.ylim([11.5,15.5])
plt.xlabel('ID [a.u.]',fontsize=16)
plt.xlim([50,600])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(numpoints=1,loc=2)
plt.grid(True)
plt.tight_layout()

plt.savefig("offset-residual-secondIter.pdf")
plt.show()
