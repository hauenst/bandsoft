#!/bin/bash
# TODO: before running, please make sure all the following paths are correct
#	(1) path2SNP: path to where the .snp files (the HV settings used in each run)
#	(2) path2Hipo: path to where the decoded .hipo files are (the actual cosmic data files)
#	(3) path2O: path to where you would like the output to go
#	(4) SNPArray: the names of the .snp files used
#	(5) HipoArray: the names of the .hipo files used
#		MAKE SURE THE ORDER MATCHES FOR (4) and (5)

# TODO: Set these appropriately:
path2SNP=/work/clas12/segarrae/bandtest/winter19-spring20/snp_settings/BAND_HV-2019_RGBSpring2019_
path2Hipo=/work/clas12/segarrae/bandtest/winter19-spring20/cosmics/hipo
path2O=/work/clas12/segarrae/bandtest/winter19-spring20/caliboutput

# TODO: Set these appropriately and make sure they match in order
SNPArray=(${path2SNP}vset_-100V_long_vset_-50V_short.snp ${path2SNP}vset_-50V_long_vset_-25V_short.snp ${path2SNP}vset_0V_long_vset_0V_short.snp ${path2SNP}vset_50V_long_vset_25V_short.snp ${path2SNP}vset_100V_long_vset_50V_short.snp)
HipoArray=(${path2Hipo}/decoded_006628.hipo ${path2Hipo}/decoded_006627.hipo ${path2Hipo}/decoded_000270.hipo ${path2Hipo}/decoded_006624.hipo ${path2Hipo}/decoded_006623.hipo)


path2OG=${path2O}/gainFitPDF
rm -r ${path2OG}
mkdir ${path2OG}
date=$(date '+%Y-%m-%d')
path2OutPNG=${path2O}/BAND_HV_${date}.snp
p2Py=../src/gainfit


for i in {0..4}
do
    #Run Spec Fit on Hipo                                                                                                                                                                    
    ./gainfit/specFit ${HipoArray[i]} ${path2O}/spectral_Fit_${i}.pdf ${path2O}/SLC_Params_${i}.txt
    #Run SNP to SLC
    python ${p2Py}/SNPtoSLC.py ${SNPArray[i]} ${path2O}/SLC_HV_${i}.txt
    #Now merge the two outputs
    ./gainfit/combFiles ${path2O}/SLC_Params_${i}.txt ${path2O}/SLC_HV_${i}.txt ${path2O}/SLC_Params_HV_${i}.txt

done

#Have a user check before continuing
echo Check the PDFs to make sure the cosmic fits are good enough.
echo After checking, manually edit the corresponding text file to have the correct ADC mean.
echo Once you complete these steps, press [Enter] to continue:
read -p Continue...

#Now fit the combined files to the fit for a gain
./gainfit/gainFit ${path2O}/HV_Settings_Order.txt ${path2OG} ${path2O}/SLC_Params_HV_*.txt

#Now restructure the output file of the fit
python ${p2Py}/combPMT.py ${path2O}/HV_Settings_Order.txt ${path2O}/HV_Settings_LR.txt

#Now make the restucted file an snp file
python ${p2Py}/mksnp.py ${path2O}/HV_Settings_LR.txt ${path2OutPNG}
