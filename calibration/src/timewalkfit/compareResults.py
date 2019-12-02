import numpy as np
import matplotlib.pyplot as plt


curr_parA_R = []; curr_parB_R = [];
prev_parA_R = []; prev_parB_R = [];
curr_parA_L = []; curr_parB_L = [];
prev_parA_L = []; prev_parB_L = [];
IDs = []

with open("../../include/time_walk_corr_right.txt","rb") as f, open("../../include/prev_time_walk_corr_right.txt","rb") as g:
	for line in f:
		parse = line.strip().split("\t")
		sector = 	int(parse[0])
		layer = 	int(parse[1])
		component = 	int(parse[2])
		barID = layer*100 + sector*10 + component
		IDs.append(barID)
		curr_parA_R.append( float(parse[3]) )
		curr_parB_R.append( float(parse[4]) )
	for line in g:
		parse = line.strip().split("\t")
		sector = 	int(parse[0])
		layer = 	int(parse[1])
		component = 	int(parse[2])
		barID = layer*100 + sector*10 + component
		prev_parA_R.append( float(parse[3]) )
		prev_parB_R.append( float(parse[4]) )


with open("../../include/time_walk_corr_left.txt","rb") as f, open("../../include/prev_time_walk_corr_left.txt","rb") as g:
	for line in f:
		parse = line.strip().split("\t")
		sector = 	int(parse[0])
		layer = 	int(parse[1])
		component = 	int(parse[2])
		barID = layer*100 + sector*10 + component
		curr_parA_L.append( float(parse[3]) )
		curr_parB_L.append( float(parse[4]) )
	for line in g:
		parse = line.strip().split("\t")
		sector = 	int(parse[0])
		layer = 	int(parse[1])
		component = 	int(parse[2])
		barID = layer*100 + sector*10 + component
		prev_parA_L.append( float(parse[3]) )
		prev_parB_L.append( float(parse[4]) )

curr_parA_R=np.asarray(curr_parA_R)
curr_parB_R=np.asarray(curr_parB_R)
curr_parA_L=np.asarray(curr_parA_L)
curr_parB_L=np.asarray(curr_parB_L)
prev_parA_R=np.asarray(prev_parA_R)
prev_parB_R=np.asarray(prev_parB_R)
prev_parA_L=np.asarray(prev_parA_L)
prev_parB_L=np.asarray(prev_parB_L)

plt.figure(1)
plt.scatter(IDs,prev_parA_R - curr_parA_R)
plt.xlim([85,665])
plt.ylim([-100,100])
plt.title('Paramter A, Right PMTs (Original-New)')
plt.xlabel('ID [a.u.]')
plt.ylabel('Par A Diff [ns]')
plt.grid(True)
plt.tight_layout()
plt.savefig('timewalk-new-right-A.pdf')

plt.figure(2)
plt.scatter(IDs,curr_parB_R - prev_parB_R)
plt.xlim([85,665])
plt.ylim([-5,5])
plt.title('Paramter B, Right PMTs (Original-New)')
plt.xlabel('ID [a.u.]')
plt.ylabel('Par A Diff [ns*sqrt(adc)]')
plt.grid(True)
plt.tight_layout()
plt.savefig('timewalk-new-right-B.pdf')

plt.figure(3)
plt.scatter(IDs,prev_parA_L - curr_parA_L)
plt.xlim([85,665])
plt.ylim([-100,100])
plt.title('Paramter A, Left PMTs (Original-New)')
plt.xlabel('ID [a.u.]')
plt.ylabel('Par A Diff [ns]')
plt.grid(True)
plt.tight_layout()
plt.savefig('timewalk-new-left-A.pdf')

plt.figure(4)
plt.scatter(IDs,curr_parB_L - prev_parB_L)
plt.xlim([85,665])
plt.ylim([-5,5])
plt.title('Paramter B, Left PMTs (Original-New)')
plt.xlabel('ID [a.u.]')
plt.ylabel('Par B Diff [ns*sqrt(adc)]')
plt.grid(True)
plt.tight_layout()
plt.savefig('timewalk-new-left-B.pdf')

#plt.show()
