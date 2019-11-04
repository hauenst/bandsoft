#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 9 09:44:45 2019

@author: adinhrnjic
"""
import decimal
import sys

if len(sys.argv)!=3:
	print 'Invalid number of arguments. Please use:\n\tpython mksnp.py [Inputfile.txt] [Outputfile.snp]'
	exit(-1)
	

infiles = [open(sys.argv[1])]

variables = [0,0,0,0,0]
# sector, layer, component, left_hv, right_hv

all_pmts = {}
layers = ["1","2","3","4","5","V"]
sides = ["_L","_R"]
numbers = []
for i in range(1,19):
    if i<10:
        numbers.append("0"+str(i))
    else:
        numbers.append(str(i))
        
default_voltage = input("Choose a default voltage for PMTs with no given setting ")

for layer in layers:
    for number in numbers:
        for side in sides:
            if int(number) > 10 and int(number) < 17:
                if layer == 'V':
                    tempa = layer + number + "A"
                    tempb = layer + number + "B"
                    all_pmts[tempa]=default_voltage
                    all_pmts[tempb]=default_voltage
                else:
                    tempa = layer + number + "A" + side
                    tempb = layer + number + "B" + side
                    all_pmts[tempa]=default_voltage
                    all_pmts[tempb]=default_voltage
            else:
                if layer == 'V':
                    temp = layer + number
                    all_pmts[temp]=default_voltage
                else:
                    temp = layer + number + side
                    all_pmts[temp]=default_voltage


for file in infiles:
    next(file)
    for line in file:
        count = 0
        for word in line.split():
            variables[count]=word
            count+=1
            
        if variables[3]=='inf':
            variables[3]=variables[4]
        if variables[4]=='inf':
            variables[4]=variables[3]
        if variables[1]=='6':
            variables[1]='V'
                        
        temp_str = 0
        if variables[0]=='1':
            if variables[1]=='V':
                temp_str = variables[1]+'0'+variables[2]
                all_pmts[temp_str]=variables[3]
            else:
                temp_str = variables[1]+'0'+variables[2]+'_L'
                all_pmts[temp_str]=variables[3]
                temp_str = variables[1]+'0'+variables[2]+'_R'
                all_pmts[temp_str]=variables[4]
                
        if variables[0]=='2':
            if variables[2]=='7':
                if variables[1]=='V':
                    temp_var = int(variables[2])+3
                    temp_str = variables[1]+str(temp_var)
                    all_pmts[temp_str]=variables[3]
                else:
                    temp_var = int(variables[2])+3
                    temp_str = variables[1]+str(temp_var)+'_L'
                    all_pmts[temp_str]=variables[3]
                    temp_str = variables[1]+str(temp_var)+'_R'
                    all_pmts[temp_str]=variables[4]
            else:
                if variables[1]=='V':
                    temp_var = int(variables[2])+3
                    temp_str = variables[1]+'0'+str(temp_var)
                    all_pmts[temp_str]=variables[3]
                else:
                    temp_var = int(variables[2])+3
                    temp_str = variables[1]+'0'+str(temp_var)+'_L'
                    all_pmts[temp_str]=variables[3]
                    temp_str = variables[1]+'0'+str(temp_var)+'_R'
                    all_pmts[temp_str]=variables[4]
                    
        if variables[0]=='3':
            if variables[1]=='V':
                temp_var = int(variables[2])+10
                temp_str = variables[1]+str(temp_var)+'A'
                all_pmts[temp_str]=variables[3]
            else:
                temp_var = int(variables[2])+10
                temp_str = variables[1]+str(temp_var)+'A_L'
                all_pmts[temp_str]=variables[3]
                temp_str = variables[1]+str(temp_var)+'A_R'
                all_pmts[temp_str]=variables[4]
                
                
        if variables[0]=='4':
            if variables[1]=='V':
                temp_var = int(variables[2])+10
                temp_str = variables[1]+str(temp_var)+'B'
                all_pmts[temp_str]=variables[4]
            else:
                temp_var = int(variables[2])+10
                temp_str = variables[1]+str(temp_var)+'B_L'
                all_pmts[temp_str]=variables[3]
                temp_str = variables[1]+str(temp_var)+'B_R'
                all_pmts[temp_str]=variables[4]

        if variables[0]=='5':
            if variables[1]=='V':
                temp_var = int(variables[2])+16
                temp_str = variables[1]+str(temp_var)
                all_pmts[temp_str]=variables[3]
            else:
                temp_var = int(variables[2])+16
                temp_str = variables[1]+str(temp_var)+'_L'
                all_pmts[temp_str]=variables[3]
                temp_str = variables[1]+str(temp_var)+'_R'
                all_pmts[temp_str]=variables[4]
            
#print(pmts)


with open(sys.argv[2],'w+') as file:
    file.write('--- Start BURT header\n')
    file.write('Time:     Tue Dec 11 17:52:57 2018\n')
    file.write('Login ID: clasrun (Online DAQ)\n')
    file.write('Eff  UID: 2508\n')
    file.write('Group ID: 9998\n')
    file.write('Keywords:\n')
    file.write('Comments:\n')
    file.write('Type:     Absolute\n')
    file.write('Directory /home/clasrun\n')
    file.write('Req File: /usr/clas12/release/1.3.0/epics/tools/burtreq/BAND_HV.req\n')
    file.write('--- End BURT header\n')
    
    for pmt in all_pmts:
        file.write('B_DET_BAND_HV_%s:vset 1 %s\n' % (pmt,'%.15e' % decimal.Decimal(all_pmts[pmt])))
        file.write('B_DET_BAND_HV_%s:vmax 1 %s\n' % (pmt,'%.15e' % decimal.Decimal(1900)))
        if( 'V' in pmt ):
		file.write('B_DET_BAND_HV_%s:iset 1 %s\n' % (pmt,'%.15e' % decimal.Decimal(500)))
	else:
		file.write('B_DET_BAND_HV_%s:iset 1 %s\n' % (pmt,'%.15e' % decimal.Decimal(300)))
        file.write('B_DET_BAND_HV_%s:trip 1 %s\n' % (pmt,'%.15e' % decimal.Decimal(0)))
        file.write('B_DET_BAND_HV_%s:rup 1 %s\n' % (pmt,'%.15e' % decimal.Decimal(50)))
        file.write('B_DET_BAND_HV_%s:rdn 1 %s\n' % (pmt,'%.15e' % decimal.Decimal(50)))
                        
