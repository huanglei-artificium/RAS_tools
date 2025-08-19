# This script cannot be run directly. It is a notebook containing the relevant commands involving the preprocessing and running the tools.

originally see /mnt/archgen/users/lei_huang/ancient_British/f4.Out.EUR_asc.2Test/*

# ascertained
for POP in Adygei Basque French Italian Orcadian Russian Sardinian Tuscan England_CNE England_WBI
do
    echo ${POP}
    Rscript ../get.blockf.R --frequencyFilesList frequency_files_list.txt --popDefFile pop_def.txt --ascOutgroup AFR_all --popLeft1 pop_list/${POP}.txt --popRight1 CEU,FIN,GBR,IBS,TSI --jackknife CHR --ascertained --ascRefpop EUR5 --ascFreqTable asc_freq_table.txt --refpopMaxMiss 0.75 --refpopRemoveHomo --outgroupMaxHetero 0 -o result/f3_block.${POP}.EUR.ascertained.txt
done

# all_sites
# need to set a larger memory size
for POP in Adygei Basque French Italian Orcadian Russian Sardinian Tuscan England_CNE England_WBI
do
    echo ${POP}
    # Rscript ../get.blockf.R --frequencyFilesList frequency_files_list.txt --ascOutgroup AFR_all --popLeft1 pop_list/${POP}.txt --popRight1 CEU,FIN,GBR,IBS,TSI --jackknife CHR -o result/f3_block.${POP}.EUR.all_sites.txt
    CMD="Rscript /mnt/archgen/users/lei_huang/RAS_tools/get.blockf.R --frequencyFilesList /mnt/archgen/users/lei_huang/RAS_tools/RAS_empirical/frequency_files_list.txt --ascOutgroup AFR_all --popLeft1 /mnt/archgen/users/lei_huang/RAS_tools/RAS_empirical/pop_list/${POP}.txt --popRight1 CEU,FIN,GBR,IBS,TSI --jackknife CHR -o /mnt/archgen/users/lei_huang/RAS_tools/RAS_empirical/result/f3_block.${POP}.EUR.all_sites.txt"
    echo ${CMD}
    qsub -V -b y -cwd -pe smp 4 -l h_vmem=128G ${CMD}
done

# 1240K
for POP in Adygei Basque French Italian Orcadian Russian Sardinian Tuscan England_CNE England_WBI
do
    echo ${POP}
    Rscript ../get.blockf.R --positionFile /mnt/archgen/users/lei_huang/genomic_data/PS.1240K.txt --frequencyFilesList frequency_files_list.txt --ascOutgroup AFR_all --popLeft1 pop_list/${POP}.txt --popRight1 CEU,FIN,GBR,IBS,TSI --jackknife CHR -o result/f3_block.${POP}.EUR.1240K.txt
done

