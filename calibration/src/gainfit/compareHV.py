from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt


if len(sys.argv)!=3:
	print 'Invalid number of arguments. Please use:\n\tpython compareHV.py [Original.snp] [New.snp]'
	exit(-1)

Map = {}
files = [sys.argv[1], sys.argv[2]]
for fi in files:
	with open(fi,'rb') as g:
		for line in g:
			if 'vset' not in line: continue
			parse = line.strip().split(" ")
			HV = float(parse[2])
			string = parse[0].split(":")[0]
			name = 0;
			if 'HV_V' in string:
				name = string.split("_")[-1]
			else:
				if 'L' in string:
					name = string.split("_")[-2]+"L"
				else:
					name = string.split("_")[-2]+"R"
					
			
			if name in Map:
				Map[name].append(HV)

			else:
				Map[name] = [HV]

IDs = []; HVs_orig = []; HVs_new = []
IDl = []; HVl_orig = []; HVl_new = []
for key in Map:
	#if 'V' not in key: continue
	if 'V' in key: 
		layer = 6
		bar = int(key[1:3])
	else:
		layer = int(key[0])
		bar = int(key[1:3])
	print  key, layer, bar

	sector = component = order = 0;

	if bar > 10 and bar < 17:
		AB = key[3]
		if layer == 6: side = ''
		else: side = key[4]

		if side == 'R': order = 1
		elif side == 'L': order = 0
		if AB == 'A': sector = 3
		elif AB == 'B': sector = 4
		
		component = bar - 10
	else:
		if layer == 6: side = ''
		else: side = key[3]

		if side == 'R': order = 1
		elif side == 'L': order = 0
		
		if bar < 4: 
			sector = 1
			component = bar
		elif bar < 11:
			sector = 2
			component = bar - 3
		elif bar < 19:
			sector = 5
			component = bar - 16
	ID = layer*1000 + sector*100 + component*10 + order
	if( len(Map[key])!=2 ): continue;
	HV_orig = Map[key][0]
	HV_new = Map[key][1]
	
	if sector == 3 or sector == 4:
		IDs.append(ID)
		HVs_orig.append( HV_orig )
		HVs_new.append( HV_new )
	else:
		IDl.append(ID)
		HVl_orig.append( HV_orig )
		HVl_new.append( HV_new )
HVs_orig = np.asarray(HVs_orig)
HVl_orig = np.asarray(HVl_orig)
HVs_new = np.asarray(HVs_new)
HVl_new = np.asarray(HVl_new)


plt.figure(1)
plt.scatter(IDs,HVs_new,color='red',label='Short Bars')
plt.scatter(IDl,HVl_new,color='blue',label='Long Bars')
plt.xlim([850,5650])
plt.ylim([-10,1600])
plt.legend(numpoints=1,loc=4)
plt.title('New HV Settings of BAND PMTs')
plt.xlabel('ID [a.u.]')
plt.ylabel('HV [V]')
plt.tight_layout()
#plt.savefig('/work/clas12/segarrae/bandtest/winter19-spring20/caliboutput/HV_new.pdf')

plt.figure(2)
plt.scatter(IDs,HVs_orig,color='red',label='Short Bars')
plt.scatter(IDl,HVl_orig,color='blue',label='Long Bars')
plt.xlim([850,5650])
plt.ylim([-10,1600])
plt.legend(numpoints=1,loc=4)
plt.title('Original HV Settings of BAND PMTs')
plt.xlabel('ID [a.u.]')
plt.ylabel('HV [V]')
plt.tight_layout()
#plt.savefig('/work/clas12/segarrae/bandtest/winter19-spring20/caliboutput/HV_old.pdf')

for i in range(len(HVs_orig)):
	if abs(HVs_orig[i] - HVs_new[i]) > 25:
		print IDs[i]
for i in range(len(HVl_orig)):
	if abs(HVl_orig[i] - HVl_new[i]) > 25:
		print IDl[i]

plt.figure(3)
plt.scatter(IDs,HVs_orig-HVs_new,color='red',label='Short Bars')
plt.scatter(IDl,HVl_orig-HVl_new,color='blue',label='Long Bars')
plt.xlim([850,6650])
plt.ylim([-100,100])
plt.legend(numpoints=1,loc=4)
plt.title('HV Differential (Original-New)')
plt.xlabel('ID [a.u.]')
plt.ylabel('HV [V]')
plt.grid(True)
plt.tight_layout()
#plt.savefig('/work/clas12/segarrae/bandtest/winter19-spring20/caliboutput/HV_diff.pdf')


plt.show()
