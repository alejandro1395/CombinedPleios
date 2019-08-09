#!/bin/bash

#set the job name
#SBATCH --job-name=Haps_Onset_Cancer
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=03:00:00

#run the application
#PATHS
OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/Combinations/results/
INPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/results/
VCF=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/data/VCF_1000G/
BIN=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/bin/plink-1.07-x86_64/
mkdir -p ${OUTPUT}CancerMetabolic_pleios/
mkdir -p ${OUTPUT}CancerMetabolic_snps/

module load Python

#### Crear haplotipos con Plink. Solo lo hace una vez por valor de R2
### Esta comanda de aqui abajo selecciona pleiotropias por el valor de r2 que tienen en la tabla.

sqlite3 ${INPUT}GWASpleiotropies.sqlite \
'SELECT DISTINCT SNPA,DiseaseA,RiskAllA,OnsetA,SNPB,DiseaseB,RiskAllB,OnsetB,ID,CHR FROM filteredPairs WHERE R2 >= 0.8 AND ((GroupA IN ("Abnormality of metabolism/homeostasis", 
 "digestive system disease", "metabolic disease") AND GroupB IN ("neoplasm", "Meningioma")) OR (GroupB IN ("Abnormality of metabolism/homeostasis", 
 "digestive system disease", "metabolic disease") AND GroupA IN ("neoplasm", "Meningioma"))) AND CHR != "" ;' > ${OUTPUT}CancerMetabolic_pleios/pleios.txt

### Construct haplotypes with plink
cat ${OUTPUT}CancerMetabolic_pleios/pleios.txt | while read line; do
    snpA=$(echo "$line" | cut -f 1);
    chr=$(echo "$line" | cut -f 10);
    echo $snpA;
    snpB=$(echo "$line" | cut -f 5);
    echo -e '*' ${snpA}'\t'${snpB} > ${OUTPUT}CancerMetabolic_snps/snps.hlist;
    sed -i -e 's/ /\t/g' ${OUTPUT}CancerMetabolic_snps/snps.hlist;
    ${BIN}plink --file ${VCF}chr${chr}_CEU_genotypes \
--hap ${OUTPUT}CancerMetabolic_snps/snps.hlist \
--hap-freq \
--noweb \
--out ${OUTPUT}CancerMetabolic_snps/${snpA}_${snpB};
grep -v LOCUS ${OUTPUT}CancerMetabolic_snps/${snpA}_${snpB}.frq.hap | awk '{print $2,$3}' > ${OUTPUT}CancerMetabolic_snps/${snpA}_${snpB}.fhtp
done

### 3 #### Create pleiotropies agon/antagon with early/late classif
cat ${OUTPUT}CancerMetabolic_pleios/pleios.txt | while read line; do
        snpA=$(echo "$line" | cut -f 1);
        snpB=$(echo "$line" | cut -f 5);
        OnsetA=$(echo "$line" | cut -f 4);
        OnsetB=$(echo "$line" | cut -f 8);
        riskHap=$(echo "$line" | awk -F"\t" '{print $3$7}');
### Python script to count pleiotropies in haplotypes
        for i in {10..60};
        do mkdir -p ${OUTPUT}CancerMetabolic_snps/Age_threeshold_${i}/
        python3 ./countPleiotropies.py ${OUTPUT}CancerMetabolic_snps/${snpA}_${snpB}.fhtp ${i} ${riskHap} ${OnsetA} ${OnsetB} ${line} ${OUTPUT}CancerMetabolic_snps/Age_threeshold_${i}/; ##### po$
    done; done

