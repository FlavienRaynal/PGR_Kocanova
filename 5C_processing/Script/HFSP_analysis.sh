 #!/bin/bash 

 # Step 1: load the requirements 

# bsub -o HFSP_analysis_29Jul14_ic_15kb1.log -e HFSP_analysis_29Jul14_ic_15kb1.err -q long -W 10:00 -R rusage[mem=400000] "bash HFSP_analysis_29Jul14_ic_15kb.sh"

# cloning the scripts 
#module load git/2.1.3
#git clone https://github.com/dekkerlab/cMapping.git
#git clone https://github.com/dekkerlab/cworld-dekker.git
#git clone https://github.com/dekkerlab/balance.git
#git clone https://github.com/blajoie/hdf2tab.git
#git clone https://github.com/blajoie/tab2hdf.git

# loading the required packages
module load python/2.7.9
module load python/2.7.9_packages/scipy/0.17.0
module load python/2.7.9_packages/numpy/1.9.2
module load python/2.7.9_packages/h5py/2.5.0
module load bedtools/2.25.0
module load R/3.3.1
## No need to "module load novocraft/V3.02.08" because you specify the path below
# alignerPath []: /home/bl73w/cMapping/alpha/aligners/novocraft/3.02.00/novoalign
# Need to create the folder below into your home directory
mkdir /home/ba69w/lsf_jobs


# Step 2: create index file
# Need to create index file for the alignment. You basically take the primer pool fasta file and convert it to a reference file for novoalign. 

# ./novoindex -m /farline/umw_job_dekker/HPCC/cshare/genome/novoalign/primer_index_name /farline/umw_job_dekker/HPCC/cshare/genome/fasta/5C/primer_index.fa
# Note : the submitted job name mapHiC should be map5C. 

 # 29JUL14_SR100_C574EAC-A  ### 6 libraries

 run_path="/farline/umw_job_dekker/HPCC/cshare/solexa/29JUL14_SR100_C574EAC-A/"
 farline_dir="/farline/umw_job_dekker/HPCC/ba69w"
 run_name="29JUL14_SR100_C574EAC-A"
 genome_index="4526-HFSPV2"
 fiveCData="/farline/umw_job_dekker/HPCC/ba69w/HFSP/fiveCData_29JUL14/filter_ic_15kb_scaled2/"
 fiveCData2="/farline/umw_job_dekker/HPCC/ba69w/HFSP/fiveCData_18Nov14/filter_ic_15kb_scaled2"
 genome_path="/farline/umw_job_dekker/HPCC/cshare/genome/"
 
  mkdir $fiveCData/singletonRemoval
  mkdir $fiveCData/singletonFiltered_anchorpurge
  mkdir $fiveCData/scaled
  mkdir $fiveCData/binmatrix
  mkdir $fiveCData/hdf2tab
  mkdir $fiveCData/balance
  mkdir $fiveCData/tab2hdf
  mkdir $fiveCData/obsMinusExp

perl ~/tools/cMapping/scripts/utilities/processFlowCell.pl \
-i $run_path \
-s $farline_dir/scratch/ \
-o $farline_dir/cData/ \
--gdir $genome_path \
-f \
--log $farline_dir/cWorld-logs/ \
--email Betul.AkgolOksuz@umassmed.edu \
-g $genome_index 


## Step 3: organize aligned files.

# removed ssh ghpcc06 from $cMapping/utilities/combine5CWrapper.sh
# If the file names are the same it rewites the file.

perl ~/tools/cMapping/scripts/utilities/combine5C.pl \
-i $farline_dir/cData/$run_name/ \
-s $farline_dir/scratch/ \
-o $farline_dir/fiveCData/ \
--gdir $genome_path \
--log $farline_dir/cWorld-logs/ \
-g $genome_index \
--short

## Step 4: Create interaction matrix 
module load bedtools/2.25.0
module load R/3.3.1

## First modify and sort the primers 

cat $farline_dir/TO_ORDER_4526_HFSPV2.txt | cut -f1,1 | grep "_FOR_" > $farline_dir/ForwardPrimers.txt
cat $farline_dir/TO_ORDER_4526_HFSPV2.txt | cut -f1,1 | grep "_REV_" > $farline_dir/ReversePrimers.txt

# sort by location
cat $farline_dir/ForwardPrimers.txt | awk -F '|' '{print $1"\t"$2"\t"$3"\t"$4}' |  awk -F ':' '{print $1"\t"$2"\t"$3"\t"$4}' | sort -k3,3 -k4,4 -k5,5 | awk '{print $1"|"$2"|"$3":"$4}' > $farline_dir/ForwardPrimers_sorted.txt 
cat $farline_dir/ReversePrimers.txt | awk -F '|' '{print $1"\t"$2"\t"$3"\t"$4}' |  awk -F ':' '{print $1"\t"$2"\t"$3"\t"$4}' | sort -k3,3 -k4,4 -k5,5 | awk '{print $1"|"$2"|"$3":"$4}' > $farline_dir/ReversePrimers_sorted.txt 

# !!! Important !!!!
# It is the first replicate but named as second. I changed the filenames from R2 to R1

# HFSP-MDAMB231-E3h-R2__4526-HFSPV2.gz
# HFSP-MCF7-E45min-R2__4526-HFSPV2.gz
# HFSP-MDAMB231-E45min-R2__4526-HFSPV2.gz
# HFSP-MCF7-E3h-R2__4526-HFSPV2.gz
# HFSP-MDAMB231-NT-R2__4526-HFSPV2.gz
# HFSP-MCF7-NT-R2__4526-HFSPV2.gz

# to 

# HFSP-MDAMB231-E3h-R1__4526-HFSPV2.gz
# HFSP-MCF7-E45min-R1__4526-HFSPV2.gz
# HFSP-MDAMB231-E45min-R1__4526-HFSPV2.gz
# HFSP-MCF7-E3h-R1__4526-HFSPV2.gz
# HFSP-MDAMB231-NT-R1__4526-HFSPV2.gz
# HFSP-MCF7-NT-R1__4526-HFSPV2.gz


# Step 5: column to matrix, replace NAs with 0 and create raw matrix heatmap

cd $run_path
for i in $(ls | grep "^HFSP") ; do \
perl ~/tools/cworld-dekker/scripts/perl/column2matrix.pl \
-i $i \
--oxh $farline_dir/ReversePrimers_sorted.txt --oyh $farline_dir/ForwardPrimers_sorted.txt -v ; done

cd $run_path
for i in $(ls | grep "score") ; do \
zcat $i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

cd $fiveCData/raw_matrix
for i in $(ls | grep "score") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $i ; done



# Step 6: Remove the diagonal. Need this step for double alternating design but not for alternating design

#cd $fiveCData/raw_matrix
#for i in $(ls | grep "score.matrix.gz") ; do \
#perl ~/tools/cworld-dekker/scripts/perl/subsetMatrix.pl \
#-i $i --minDist 1 -v ; done

#cd $fiveCData/raw_matrix
#for i in $(ls | grep "subset") ; do \
#perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
#-i $i ; done


# Step 7: remove singletons 

input=$fiveCData/raw_matrix
output=$fiveCData/singletonRemoval
cd $output
for i in $(ls $input | grep "score.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/singletonRemoval.pl \
-i $input/$i --ic --it --ca 0.02 --caf 2500 --tta 6 --cta 12 --ez  -v ; done

# Optional 
# remove singletions again using the union of removed singletons for all replicates 

# take the union of singletons
cat $fiveCData/singletonRemoval/*.toRemove.txt \
$fiveCData/singletonRemoval/*.toRemove.txt \
| awk '{print $2}'| awk '!seen[$0]++' > $fiveCData/HFSP_all_samples.singletonRemoval.toRemove.txt

# remove union singletons 
input=$fiveCData/raw_matrix
output=$fiveCData/singletonRemoval
cd $output
for i in $(ls $input | grep "score.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/singletonRemoval.pl \
-i $input/$i --ic --it --ca 0.02 --caf 2500 --tta 6 --cta 12 --ez  --mof $fiveCData/HFSP_all_samples.singletonRemoval.toRemove.txt -v ; done


# Step 8: Remove anchors

input=$fiveCData/singletonRemoval/
output=$fiveCData/singletonFiltered_anchorpurge
cd $output
for i in $(ls $input | grep "singletonFiltered.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/anchorPurge.pl \
-i $input/$i --ic --ca 0.02 --caf 2500 -v ; done

# take the union of anchorpurge 
cat $fiveCData/singletonFiltered_anchorpurge/*.toRemove \
$fiveCData2/singletonFiltered_anchorpurge/*.toRemove \
| awk '{print $1}'| awk '!seen[$0]++' > $fiveCData/HFSP_all_samples_anchorpurge.toRemove.txt

# remove union anchorpurge 
input=$fiveCData/singletonRemoval/
output=$fiveCData/singletonFiltered_anchorpurge
cd $output
for i in $(ls $input | grep "singletonFiltered.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/anchorPurge.pl \
-i $input/$i --ic --ca 0.02 --caf 2500 --mof $fiveCData/HFSP_all_anchorpurge.toRemove.txt -v ; done



# Step 9: Make the matrix symmetrical

input=$fiveCData/singletonFiltered_anchorpurge/
output=$fiveCData/singletonFiltered_anchorpurge
cd $output
for i in $(ls | grep "singletonFiltered.anchorFiltered.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/matrix2symmetrical.pl \
-i $i -v ; done


# Step 10: scale the matrix to 15000000

input=$fiveCData/singletonFiltered_anchorpurge/
output=$fiveCData/scaled
cd $output
for i in $(ls $input | grep "singletonFiltered.anchorFiltered.symmetrical.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/scaleMatrix.pl \
-i $input/$i --st 15000000 -v ; done

# # Step 11: Balancing 1 
input=$fiveCData/scaled/
output=$fiveCData/tab2hdf/
cd $output
for i in $(ls $input | grep "15000000.matrix.gz") ; do \
python ~/tools/tab2hdf/scripts/tab2hdf.py \
-i $input/$i \
-v ; done


input=$fiveCData/tab2hdf/
output=$fiveCData/balance/
cd $output
for i in $(ls $input | grep "symmetrical.scaled-15000000.hdf5") ; do \
python ~/tools/balance/scripts/balance.py  \
-i $input/$i \
-v ; done

input=$fiveCData/balance/
output=$fiveCData/hdf2tab/
cd $output
for i in $(ls $input | grep "balanced.hdf5") ; do \
python ~/tools/hdf2tab/scripts/hdf2tab.py \
-i $input/$i \
-v ; done



# Step 11: Binning

input=$fiveCData/hdf2tab/
output=$fiveCData/binmatrix/
cd $output
for i in $(ls $input | grep ".balanced.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/binMatrix.pl --bsize 15000 --bstep 8 --bmode median \
-i $input/$i ; done


#  Step 12: Balancing 2

input=$fiveCData/binmatrix/
output=$fiveCData/tab2hdf/
cd $output
for i in $(ls $input | grep "__8.matrix.gz") ; do \
python ~/tools/tab2hdf/scripts/tab2hdf.py \
-i $input/$i \
-v ; done

input=$fiveCData/tab2hdf/
output=$fiveCData/balance/
cd $output
for i in $(ls $input | grep "__8.hdf5") ; do \
python ~/tools/balance/scripts/balance.py  \
-i $input/$i \
-v ; done

input=$fiveCData/balance/
output=$fiveCData/hdf2tab/
cd $output
for i in $(ls $input | grep "__8.balanced.hdf5") ; do \
python ~/tools/hdf2tab/scripts/hdf2tab.py \
-i $input/$i \
-v ; done

# Step 13: Observed minus expected 

input=$fiveCData/hdf2tab/
output=$fiveCData/obsMinusExp/
cd $output
for i in $(ls $input | grep "balanced.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/matrix2loess.pl \
-i $input/$i \
-v ; done


# Step 14: Compare treatmen to non treatment and two cell lines 
# 1. compare NT to 3h for both cell lines
# 2. Compare MCF7 NT to MDA-MD231 NT

# Using balance matrix
 MCF7_NT="HFSP-MCF7-NT-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.matrix.gz"
 MDA_MD231_NT="HFSP-MDAMB231-NT-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.matrix.gz"
 MCF7_3h="HFSP-MCF7-E3h-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.matrix.gz"
 MDA_MD231_3h="HFSP-MDAMB231-E3h-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.matrix.gz"

 input=$fiveCData/hdf2tab/
 output=$fiveCData/comparison
 cd $output
 perl ~/tools/cworld-dekker/scripts/perl/compareMatrices.pl \
 -1 $input/$MCF7_NT -2 $input/$MCF7_3h --cm log2ratio -o MCF7_NT_MCF7_3h

 perl ~/tools/cworld-dekker/scripts/perl/compareMatrices.pl \
 -1 $input/$MDA_MD231_NT -2 $input/$MDA_MD231_3h --cm log2ratio -o MDA_MD231_NT_MDA_MD231_3h

 perl ~/tools/cworld-dekker/scripts/perl/compareMatrices.pl \
 -1 $input/$MCF7_NT -2 $input/$MDA_MD231_NT --cm log2ratio -o MCF7_NT_MDA_MD231_NT


#  Heatmap of step 14
 perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
 -i $output/MCF7_NT_MCF7_3h.log2ratio.matrix.gz

 perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
 -i $output/MDA_MD231_NT_MDA_MD231_3h.log2ratio.matrix.gz 

 perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
 -i $output/MCF7_NT_MDA_MD231_NT.log2ratio.matrix.gz


# Using obs/exp matrix
 MDA_MD231_NT="HFSP-MDAMB231-NT-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.obs-exp.matrix.gz"
 MDA_MD231_3h="HFSP-MDAMB231-E3h-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.obs-exp.matrix.gz"
 MCF7_NT="HFSP-MCF7-NT-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.obs-exp.matrix.gz"
 MCF7_3h="HFSP-MCF7-E3h-R1__4526-HFSPV2.score.singletonFiltered.anchorFiltered.symmetrical.scaled-15000000.balanced__15000__8.balanced.obs-exp.matrix.gz"

 input=$fiveCData/obsMinusExp/
 output=$fiveCData/comparison2
 cd $output
 perl ~/tools/cworld-dekker/scripts/perl/compareMatrices.pl \
 -1 $input/$MCF7_NT -2 $input/$MCF7_3h --cm subtract -o MCF7_NT_MCF7_3h

 perl ~/tools/cworld-dekker/scripts/perl/compareMatrices.pl \
 -1 $input/$MDA_MD231_NT -2 $input/$MDA_MD231_3h --cm subtract -o MDA_MD231_NT_MDA_MD231_3h

 perl ~/tools/cworld-dekker/scripts/perl/compareMatrices.pl \
 -1 $input/$MCF7_NT -2 $input/$MDA_MD231_NT --cm subtract -o MCF7_NT_MDA_MD231_NT

# # Heatmap of obs/exp matrix
 perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
 -i $output/MCF7_NT_MCF7_3h.subtract.matrix.gz

 perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
 -i $output/MDA_MD231_NT_MDA_MD231_3h.subtract.matrix.gz 

 perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
 -i $output/MCF7_NT_MDA_MD231_NT.subtract.matrix.gz



###########  Heatmaps of all steps    ###################


# Heatmap singleton removed matrix
input=$fiveCData/singletonRemoval/
output=$fiveCData/singletonRemoval
cd $output
for i in $(ls $input | grep "singletonFiltered.matrix.gz") ; do \
zcat $input/$i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

input=$fiveCData/singletonRemoval/
output=$fiveCData/singletonRemoval
cd $output
for i in $(ls $input | grep "singletonFiltered.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $input/$i ; done

# Remove anchors 
input=$fiveCData/singletonFiltered_anchorpurge/
output=$fiveCData/singletonFiltered_anchorpurge
cd $output
for i in $(ls $input | grep "singletonFiltered.anchorFiltered.matrix.gz") ; do \
zcat $input/$i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

input=$fiveCData/singletonFiltered_anchorpurge/
output=$fiveCData/singletonFiltered_anchorpurge
cd $output
for i in $(ls $input | grep "singletonFiltered.anchorFiltered.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $input/$i ; done


# Heatmap make symmetrical
input=$fiveCData/singletonFiltered_anchorpurge/
output=$fiveCData/singletonFiltered_anchorpurge
cd $output
for i in $(ls $input | grep "singletonFiltered.anchorFiltered.symmetrical.matrix.gz") ; do \
zcat $input/$i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

input=$fiveCData/singletonFiltered_anchorpurge/
output=$fiveCData/singletonFiltered_anchorpurge
cd $output
for i in $(ls $input | grep "singletonFiltered.anchorFiltered.symmetrical.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $input/$i ; done

# HEATMAP OF scaled matrix 

input=$fiveCData/scaled/
output=$fiveCData/scaled
cd $output
for i in $(ls $input | grep ".matrix.gz") ; do \
zcat $input/$i| awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

input=$fiveCData/scaled/
output=$fiveCData/scaled
cd $output
for i in $(ls $input | grep ".matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $input/$i ; done


# HEATMAP OF binned matrix 

input=$fiveCData/binmatrix/
output=$fiveCData/binmatrix
cd $output
for i in $(ls $input | grep ".matrix.gz") ; do \
zcat $input/$i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

input=$fiveCData/binmatrix/
output=$fiveCData/binmatrix
cd $output
for i in $(ls $input | grep ".matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $input/$i ; done



# Balace 1
# Balancing 2
# Replace "nan" with 0

#input=$fiveCData/hdf2tab/
#output=$fiveCData/hdf2tab
#for i in $(ls | grep "symmetrical.balanced.matrix.gz") ; do \
#zcat $i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done


input=$fiveCData/hdf2tab/
output=$fiveCData/hdf2tab
cd $output
for i in $(ls $input | grep "matrix.gz") ; do \
zcat $input/$i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^nan*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

input=$fiveCData/hdf2tab/
output=$fiveCData/hdf2tab
cd $output
for i in $(ls $input | grep ".matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $input/$i ; done


# Heatmap observed/expected

input=$fiveCData/obsMinusExp/
output=$fiveCData/obsMinusExp
cd $output
for i in $(ls $input | grep "obs-exp.matrix.gz") ; do \
zcat $input/$i | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^NA*$/) $i = 0 }; 1' | gzip > "$i.mod.txt.gz" ; done

input=$fiveCData/obsMinusExp/
output=$fiveCData/obsMinusExp
cd $output
for i in $(ls $input | grep ".obs-exp.matrix.gz") ; do \
perl ~/tools/cworld-dekker/scripts/perl/heatmap.pl \
-i $input/$i ; done





