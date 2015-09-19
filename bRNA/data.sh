#!/bin/sh
# AUTHOR: XULONG WANG (XULONG.WANG@JAX.ORG)
# FUNCTION: to make data matrix from cuffdiff/rsem outputs
# USEAGE: mydata.sh

## Cuffdiff output
#cd ./cuffdiff_mm10
#
#for myfile in genes.read_group_tracking isoforms.read_group_tracking 
#do
#awk '{if ($2=="q1" && $3=="0") print $1,$4}' $myfile > vninp1
#awk '{if ($2=="q1" && $3=="1") print $1,$4}' $myfile > vnipp1
#awk '{if ($2=="q1" && $3=="2") print $1,$4}' $myfile > vpipp1
#awk '{if ($2=="q2" && $3=="0") print $1,$4}' $myfile > vninp2
#awk '{if ($2=="q2" && $3=="1") print $1,$4}' $myfile > vnipp2
#awk '{if ($2=="q2" && $3=="2") print $1,$4}' $myfile > vpipp2
#awk '{if ($2=="q1" && $3=="1") print $1,$4}' $inputfile > q1_1
#awk '{if ($2=="q1" && $3=="2") print $1,$4}' $inputfile > q1_2
#awk '{if ($2=="q1" && $3=="3") print $1,$4}' $inputfile > q1_3
#
#awk '{if ($2=="q2" && $3=="0") print $1,$4}' $inputfile > q2_0
#awk '{if ($2=="q2" && $3=="1") print $1,$4}' $inputfile > q2_1
#awk '{if ($2=="q2" && $3=="2") print $1,$4}' $inputfile > q2_2
#awk '{if ($2=="q2" && $3=="3") print $1,$4}' $inputfile > q2_3
#awk '{if ($2=="q2" && $3=="4") print $1,$4}' $inputfile > q2_4
#
#join vninp1 vninp2 > t1
#join vnipp1 vnipp2 > t2
#join vpipp1 vpipp2 > t3
#join t1 t2 > t4
#join t4 t3 > ${myfile:-s/read_group_tracking/output}
#join t4 t3 > "$myfile"_output
#done
#
#join q2_0 q2_1 > t1
#join t1 q2_2 > t2
#join t2 q2_3 > t3
#join t3 q2_4 > q2
#
##join q1 q2 > output.genes
#join q1 q2 > output.isoforms
#
#rm t[123]
#rm q*

### RSEM output

cd ~/Dropbox/Lupus/RSEM
#for myfile in *.results
#do
#   awk '{print $1,$5}' $myfile > "$myfile"_count
#	awk '{print $1,$6}' $myfile > "$myfile"_TPM
#	awk '{print $1,$7}' $myfile > "$myfile"_FPKM
#done
#
type="TPM"  # count, TPM, FPKM
for myfile in genes isoforms
do
    join VnInp1."$myfile".results_"$type" VnInp2."$myfile".results_"$type" > t1
    join t1 VnIpp1."$myfile".results_"$type" > t2
    join t2 VnIpp2."$myfile".results_"$type" > t3
    join t3 VpIpp1."$myfile".results_"$type" > t4
    join t4 VpIpp2."$myfile".results_"$type" > "$myfile"."$type"
done
#
rm t[1234]

myinput=$1
# Input file: one column of Ensembl IDs 

join $myinput ensembl2genesymbol.txt > temp
awk '{print $2}' temp > $myinput

