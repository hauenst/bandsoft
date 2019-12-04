from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

s_veff_tdc = [] ; 
s_veff_ftdc = []; 
s_off_tdc = []; 
s_off_ftdc = [];
l_veff_tdc = [] ; 
l_veff_ftdc = []; 
l_off_tdc = []; 
l_off_ftdc = [];
IDs = [];
veff_tdc = [];
veff_ftdc = [];
off_tdc = [];
off_ftdc = [];
with open("../../include/effective_velocity.txt","rb") as f, open("../../include/lr_offsets.txt","rb") as g:
	for line in f:
		parse = line.strip().split("\t")
		if( float(parse[3])==0 or float(parse[4])==0): continue
		ID = int(parse[1])*100 + int(parse[0])*10 + int(parse[2])
		IDs.append(ID); veff_tdc.append( float(parse[3]) ); veff_ftdc.append( float(parse[4]) )
		if( int(parse[0]) == 3 or int(parse[0]) == 4):
			s_veff_tdc.append( float(parse[3]) )
			s_veff_ftdc.append( float(parse[4]) )
		else:
			l_veff_tdc.append( float(parse[3]) )
			l_veff_ftdc.append( float(parse[4]) )
	for line in g:
		parse = line.strip().split("\t")
		if( parse[3]==0 or parse[4]==0): continue
		ID = int(parse[1])*100 + int(parse[0])*10 + int(parse[2])
		off_tdc.append( float(parse[3]) ); off_ftdc.append( float(parse[4]) )
		if( int(parse[0]) == 3 or int(parse[0]) == 4):
			s_off_tdc.append( float(parse[3]) )
			s_off_ftdc.append( float(parse[4]) )
		else:
			l_off_tdc.append( float(parse[3]) )
			l_off_ftdc.append( float(parse[4]) )

plt.figure(1)
plt.scatter(s_veff_tdc,s_veff_ftdc,color='red',label='Short bars')
plt.scatter(l_veff_tdc,l_veff_ftdc,color='blue',label='Long bars')
plt.xlabel("TDC-Measured Effective Velocity [cm/ns]",fontsize=16)
plt.ylabel("FADC-Measured Effective Velocity [cm/ns]",fontsize=16)
plt.xlim([11.5,15.5])
plt.ylim([11.5,15.5])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(numpoints=1,loc=2)
plt.grid(True)
plt.tight_layout()
plt.savefig("eff_vel_tdcfadc.pdf")

plt.figure(2)
plt.scatter(IDs,veff_tdc,color='red',label='TDC')
plt.scatter(IDs,veff_ftdc,color='blue',label='FADC')
plt.ylabel("Effective Velocity [cm/ns]",fontsize=16)
plt.ylim([11.5,15.5])
plt.xlabel('ID [a.u.]',fontsize=16)
plt.xlim([50,600])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(numpoints=1,loc=2)
plt.grid(True)
plt.tight_layout()
plt.savefig("eff_vel_individual.pdf")


plt.show()
