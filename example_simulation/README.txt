# original data
/mnt/archgen/MICROSCOPE/rarevar_sim_study/data/from_mpcdf/n50/10000000

cd /mnt/archgen/users/lei_huang/RAS_tools/example_simulation

mkdir data
for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    mkdir data/n50_${RATE}
done


cd trash
# will generate lots of cluster files. run in trash folder
# because of dump or malfunction on some nodes, some tasks will fail to finish. need to run multiple times (judge if the file already exists before generating a file)
for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    echo ${RATE}
    original_dir="/mnt/archgen/MICROSCOPE/rarevar_sim_study/data/from_mpcdf/n50/10000000/${RATE}.0"
    output_dir="/mnt/archgen/users/lei_huang/RAS_tools/example_simulation/data/n50_${RATE}"

    for POP_IDX in {0..8}
    do
        OFFSET=$[50*POP_IDX]

        for GRP_IDX in {0..4}
        do
            echo ${POP_IDX} ${GRP_IDX}
            GRP_ST=$[OFFSET+10*GRP_IDX+1]
            GRP_ED=$[OFFSET+10*GRP_IDX+10]
            echo ${GRP_ST} ${GRP_ED}

            if [ ! -f ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.geno ]
            then
                # cut -c ${GRP_ST}-${GRP_ED} ${original_dir}/all_vars.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.geno
                CMD="cut -c ${GRP_ST}-${GRP_ED} ${original_dir}/all_vars.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.geno"
                echo ${CMD}
                qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
            fi
            if [ ! -f ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.geno ]
            then
                # cut -c ${GRP_ST}-${GRP_ED} ${original_dir}/twelve_forty_vars.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.geno
                CMD="cut -c ${GRP_ST}-${GRP_ED} ${original_dir}/twelve_forty_vars.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.geno"
                echo ${CMD}
                qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
            fi
            if [ ! -f ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.ind ]
            then
                # sed -n "${GRP_ST},${GRP_ED}p" ${original_dir}/all_vars.ind > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.ind
                CMD="sed -n \"${GRP_ST},${GRP_ED}p\" ${original_dir}/all_vars.ind > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.ind"
                echo ${CMD}
                qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
            fi
        done
    done

    if [ ! -f ${output_dir}/out.all_vars.geno ]
    then
        # cut -c 451 ${original_dir}/all_vars.geno > ${output_dir}/out.all_vars.geno
        CMD="cut -c 451 ${original_dir}/all_vars.geno > ${output_dir}/out.all_vars.geno"
        echo ${CMD}
        qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
    fi
    if [ ! -f ${output_dir}/out.1240K.geno ]
    then
        # cut -c 451 ${original_dir}/twelve_forty_vars.geno > ${output_dir}/out.1240K.geno
        CMD="cut -c 451 ${original_dir}/twelve_forty_vars.geno > ${output_dir}/out.1240K.geno"
        echo ${CMD}
        qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
    fi
    if [ ! -f ${output_dir}/out.ind ]
    then
        # sed -n "451p" ${original_dir}/all_vars.ind > ${output_dir}/out.ind
        CMD="sed -n \"451p\" ${original_dir}/all_vars.ind > ${output_dir}/out.ind"
        echo ${CMD}
        qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
    fi

    if [ ! -f ${output_dir}/PS.all_vars.txt ]
    then
        # cut -f 2,4 ${original_dir}/all_vars.snp > ${output_dir}/PS.all_vars.txt
        CMD="cut -f 2,4 ${original_dir}/all_vars.snp > ${output_dir}/PS.all_vars.txt"
        echo ${CMD}
        qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
    fi
    if [ ! -f ${output_dir}/PS.1240K.txt ]
    then
        # cut -f 2,4 ${original_dir}/twelve_forty_vars.snp > ${output_dir}/PS.1240K.txt
        CMD="cut -f 2,4 ${original_dir}/twelve_forty_vars.snp > ${output_dir}/PS.1240K.txt"
        echo ${CMD}
        qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
    fi

done


cd trash
for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    echo ${RATE}
    output_dir="/mnt/archgen/users/lei_huang/RAS_tools/example_simulation/data/n50_${RATE}"

    for POP_IDX in {0..8}
    do
        for GRP_IDX in {0..4}
        do
            echo ${POP_IDX} ${GRP_IDX}

            if [ ! -f ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.frac.txt ]
            then
                # awk -F '' -v OFS='\t' '{n=0; sum=0; for (i=1;i<=NF;i++) if ($i!=9) {n++; sum+=$i} print sum,n*2}' ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.frac.txt
                CMD="awk -F '' -v OFS='\t' '{n=0; sum=0; for (i=1;i<=NF;i++) if (\$i!=9) {n++; sum+=\$i} print sum,n*2}' ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.all_vars.frac.txt"
                echo ${CMD}
                qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
            fi
            if [ ! -f ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.frac.txt ]
            then
                # awk -F '' -v OFS='\t' '{n=0; sum=0; for (i=1;i<=NF;i++) if ($i!=9) {n++; sum+=$i} print sum,n*2}' ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.frac.txt
                CMD="awk -F '' -v OFS='\t' '{n=0; sum=0; for (i=1;i<=NF;i++) if (\$i!=9) {n++; sum+=\$i} print sum,n*2}' ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.geno > ${output_dir}/pop${POP_IDX}.group${GRP_IDX}.1240K.frac.txt"
                echo ${CMD}
                qsub -V -b y -cwd -pe smp 4 -l h_vmem=16G ${CMD}
            fi
        done
    done
done



for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    echo ${RATE}
    output_dir="/mnt/archgen/users/lei_huang/RAS_tools/example_simulation/data/n50_${RATE}"

    for POP_IDX in {0..8}
    do
        echo ${POP_IDX}
        if [ ! -f ${output_dir}/pop${POP_IDX}.x_groupi.all_vars.frac.txt ]
        then
            paste ${output_dir}/pop${POP_IDX}.group0.all_vars.frac.txt ${output_dir}/pop${POP_IDX}.group1.all_vars.frac.txt ${output_dir}/pop${POP_IDX}.group2.all_vars.frac.txt ${output_dir}/pop${POP_IDX}.group3.all_vars.frac.txt ${output_dir}/pop${POP_IDX}.group4.all_vars.frac.txt | awk -F '\t' -v OFS='\t' '{REF_ALL=$1+$3+$5+$7+$9; COUNT_ALL=$2+$4+$6+$8+$10; print REF_ALL-$1,COUNT_ALL-$2,REF_ALL-$3,COUNT_ALL-$4,REF_ALL-$5,COUNT_ALL-$6,REF_ALL-$7,COUNT_ALL-$8,REF_ALL-$9,COUNT_ALL-$10}' > ${output_dir}/pop${POP_IDX}.x_groupi.all_vars.frac.txt
        fi
        if [ ! -f ${output_dir}/pop${POP_IDX}.x_groupi.1240K.frac.txt ]
        then
            paste ${output_dir}/pop${POP_IDX}.group0.1240K.frac.txt ${output_dir}/pop${POP_IDX}.group1.1240K.frac.txt ${output_dir}/pop${POP_IDX}.group2.1240K.frac.txt ${output_dir}/pop${POP_IDX}.group3.1240K.frac.txt ${output_dir}/pop${POP_IDX}.group4.1240K.frac.txt | awk -F '\t' -v OFS='\t' '{REF_ALL=$1+$3+$5+$7+$9; COUNT_ALL=$2+$4+$6+$8+$10; print REF_ALL-$1,COUNT_ALL-$2,REF_ALL-$3,COUNT_ALL-$4,REF_ALL-$5,COUNT_ALL-$6,REF_ALL-$7,COUNT_ALL-$8,REF_ALL-$9,COUNT_ALL-$10}' > ${output_dir}/pop${POP_IDX}.x_groupi.1240K.frac.txt
        fi
    done
done


for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    echo ${RATE}
    output_dir="/mnt/archgen/users/lei_huang/RAS_tools/example_simulation/data/n50_${RATE}"

    paste ${output_dir}/pop0.x_groupi.all_vars.frac.txt ${output_dir}/pop1.x_groupi.all_vars.frac.txt ${output_dir}/pop2.x_groupi.all_vars.frac.txt ${output_dir}/pop3.x_groupi.all_vars.frac.txt ${output_dir}/pop4.x_groupi.all_vars.frac.txt ${output_dir}/pop5.x_groupi.all_vars.frac.txt ${output_dir}/pop6.x_groupi.all_vars.frac.txt ${output_dir}/pop7.x_groupi.all_vars.frac.txt ${output_dir}/pop8.x_groupi.all_vars.frac.txt | awk -F '\t' -v OFS='\t' '{print $1+$11+$21+$31+$41+$51+$61+$71+$81,$2+$12+$22+$32+$42+$52+$62+$72+$82,$3+$13+$23+$33+$43+$53+$63+$73+$83,$4+$14+$24+$34+$44+$54+$64+$74+$84,$5+$15+$25+$35+$45+$55+$65+$75+$85,$6+$16+$26+$36+$46+$56+$66+$76+$86,$7+$17+$27+$37+$47+$57+$67+$77+$87,$8+$18+$28+$38+$48+$58+$68+$78+$88,$9+$19+$29+$39+$49+$59+$69+$79+$89,$10+$20+$30+$40+$50+$60+$70+$80+$90}' > ${output_dir}/pop_all.x_groupi.all_vars.frac.txt

    paste ${output_dir}/pop0.x_groupi.1240K.frac.txt ${output_dir}/pop1.x_groupi.1240K.frac.txt ${output_dir}/pop2.x_groupi.1240K.frac.txt ${output_dir}/pop3.x_groupi.1240K.frac.txt ${output_dir}/pop4.x_groupi.1240K.frac.txt ${output_dir}/pop5.x_groupi.1240K.frac.txt ${output_dir}/pop6.x_groupi.1240K.frac.txt ${output_dir}/pop7.x_groupi.1240K.frac.txt ${output_dir}/pop8.x_groupi.1240K.frac.txt | awk -F '\t' -v OFS='\t' '{print $1+$11+$21+$31+$41+$51+$61+$71+$81,$2+$12+$22+$32+$42+$52+$62+$72+$82,$3+$13+$23+$33+$43+$53+$63+$73+$83,$4+$14+$24+$34+$44+$54+$64+$74+$84,$5+$15+$25+$35+$45+$55+$65+$75+$85,$6+$16+$26+$36+$46+$56+$66+$76+$86,$7+$17+$27+$37+$47+$57+$67+$77+$87,$8+$18+$28+$38+$48+$58+$68+$78+$88,$9+$19+$29+$39+$49+$59+$69+$79+$89,$10+$20+$30+$40+$50+$60+$70+$80+$90}' > ${output_dir}/pop_all.x_groupi.1240K.frac.txt
done


----------------------

Rscript ../get.blockf.R --frequencyFilesList frequency_files_list/n50_1.1240K.txt --popLeft1 pop0.group4,pop1.group4,pop2.group4,pop3.group4,pop4.group4,pop5.group4,pop6.group4,pop7.group4,pop8.group4 --popRight1 pop0.x_group4,pop1.x_group4,pop2.x_group4,pop3.x_group4,pop4.x_group4,pop5.x_group4,pop6.x_group4,pop7.x_group4,pop8.x_group4 --jackknife CHR --ascOutgroup Ref --allSitesBlock -o result/n50_1/f3_block.group4.x_group4.1240K.txt



for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    for GRP_IDX in 0 1 2 3 4
    do
        Rscript ../get.blockf.R --frequencyFilesList frequency_files_list/n50_${RATE}.1240K.txt --popLeft1 pop0.group${GRP_IDX},pop1.group${GRP_IDX},pop2.group${GRP_IDX},pop3.group${GRP_IDX},pop4.group${GRP_IDX},pop5.group${GRP_IDX},pop6.group${GRP_IDX},pop7.group${GRP_IDX},pop8.group${GRP_IDX} --popRight1 pop0.x_group${GRP_IDX},pop1.x_group${GRP_IDX},pop2.x_group${GRP_IDX},pop3.x_group${GRP_IDX},pop4.x_group${GRP_IDX},pop5.x_group${GRP_IDX},pop6.x_group${GRP_IDX},pop7.x_group${GRP_IDX},pop8.x_group${GRP_IDX} --jackknife CHR --ascOutgroup Ref --allSitesBlock -o result/n50_${RATE}/f3_block.group${GRP_IDX}.x_group${GRP_IDX}.1240K.txt
    done
done


for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    for GRP_IDX in 0 1 2 3 4
    do
        Rscript ../get.blockf.R --frequencyFilesList frequency_files_list/n50_${RATE}.all_vars.txt --popLeft1 pop0.group${GRP_IDX},pop1.group${GRP_IDX},pop2.group${GRP_IDX},pop3.group${GRP_IDX},pop4.group${GRP_IDX},pop5.group${GRP_IDX},pop6.group${GRP_IDX},pop7.group${GRP_IDX},pop8.group${GRP_IDX} --popRight1 pop0.x_group${GRP_IDX},pop1.x_group${GRP_IDX},pop2.x_group${GRP_IDX},pop3.x_group${GRP_IDX},pop4.x_group${GRP_IDX},pop5.x_group${GRP_IDX},pop6.x_group${GRP_IDX},pop7.x_group${GRP_IDX},pop8.x_group${GRP_IDX} --jackknife CHR --ascOutgroup Ref --allSitesBlock -o result/n50_${RATE}/f3_block.group${GRP_IDX}.x_group${GRP_IDX}.all_vars.txt
    done
done


for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    for GRP_IDX in 0 1 2 3 4
    do
        Rscript ../get.blockf.R --frequencyFilesList frequency_files_list/n50_${RATE}.all_vars.txt --popLeft1 pop0.group${GRP_IDX},pop1.group${GRP_IDX},pop2.group${GRP_IDX},pop3.group${GRP_IDX},pop4.group${GRP_IDX},pop5.group${GRP_IDX},pop6.group${GRP_IDX},pop7.group${GRP_IDX},pop8.group${GRP_IDX} --popRight1 pop0.x_group${GRP_IDX},pop1.x_group${GRP_IDX},pop2.x_group${GRP_IDX},pop3.x_group${GRP_IDX},pop4.x_group${GRP_IDX},pop5.x_group${GRP_IDX},pop6.x_group${GRP_IDX},pop7.x_group${GRP_IDX},pop8.x_group${GRP_IDX} --jackknife CHR --ascertained --ascRefpop pop_all.x_group${GRP_IDX} --ascFreqTable asc_freq_table.txt --allSitesBlock --refpopMaxMiss 0.75 --refpopRemoveHomo --ascOutgroup Ref --outgroupMaxHetero 0 -o result/n50_${RATE}/f3_block.group${GRP_IDX}.x_group${GRP_IDX}.ascertained.txt
    done
done


# if run separately, then there could be such error:
# libgomp: Thread creation failed: Resource temporarily unavailable
# qsub -V -b y -cwd -pe smp 8 -l h_vmem=128G

for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    for POP_IDX in 0 1 2 3 4 5 6 7 8
    do
        for GRP_IDX in 0 1 2 3 4
        do
            if [ ! -f /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/result/n50_${RATE}/f3_block.group${GRP_IDX}_pop${POP_IDX}_ind.x_group${GRP_IDX}.ascertained.txt ]
            then
                Rscript /mnt/archgen/users/lei_huang/RAS_tools/get.blockf.R --frequencyFilesList /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/frequency_files_list/n50_${RATE}.all_vars.txt --popLeft1 <(cut -f 1 /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/data/n50_${RATE}/pop${POP_IDX}.group${GRP_IDX}.ind) --popRight1 pop0.x_group${GRP_IDX},pop1.x_group${GRP_IDX},pop2.x_group${GRP_IDX},pop3.x_group${GRP_IDX},pop4.x_group${GRP_IDX},pop5.x_group${GRP_IDX},pop6.x_group${GRP_IDX},pop7.x_group${GRP_IDX},pop8.x_group${GRP_IDX} --jackknife CHR --ascertained --ascRefpop pop_all.x_group${GRP_IDX} --ascFreqTable /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/asc_freq_table.txt --allSitesBlock --refpopMaxMiss 0.75 --refpopRemoveHomo --ascOutgroup Ref --outgroupMaxHetero 0 -o /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/result/n50_${RATE}/f3_block.group${GRP_IDX}_pop${POP_IDX}_ind.x_group${GRP_IDX}.ascertained.txt
            fi
        done
    done
done




for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    for POP_IDX in 0 1 2 3 4 5 6 7 8
    do
        for GRP_IDX in 0 1 2 3 4
        do
            if [ ! -f /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/result/n50_${RATE}/f3_block.group${GRP_IDX}_pop${POP_IDX}_ind.x_group${GRP_IDX}.all_vars.txt ]
            then
                Rscript /mnt/archgen/users/lei_huang/RAS_tools/get.blockf.R --frequencyFilesList /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/frequency_files_list/n50_${RATE}.all_vars.txt --popLeft1 <(cut -f 1 /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/data/n50_${RATE}/pop${POP_IDX}.group${GRP_IDX}.ind) --popRight1 pop0.x_group${GRP_IDX},pop1.x_group${GRP_IDX},pop2.x_group${GRP_IDX},pop3.x_group${GRP_IDX},pop4.x_group${GRP_IDX},pop5.x_group${GRP_IDX},pop6.x_group${GRP_IDX},pop7.x_group${GRP_IDX},pop8.x_group${GRP_IDX} --jackknife CHR --allSitesBlock --ascOutgroup Ref --outgroupMaxHetero 0 -o /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/result/n50_${RATE}/f3_block.group${GRP_IDX}_pop${POP_IDX}_ind.x_group${GRP_IDX}.all_vars.txt
            fi
        done
    done
done




for RATE in 1 2 5 10 20 50 100 200 500 1000 2000 5000
do
    for POP_IDX in 0 1 2 3 4 5 6 7 8
    do
        for GRP_IDX in 0 1 2 3 4
        do
            if [ ! -f /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/result/n50_${RATE}/f3_block.group${GRP_IDX}_pop${POP_IDX}_ind.x_group${GRP_IDX}.1240K.txt ]
            then
                Rscript /mnt/archgen/users/lei_huang/RAS_tools/get.blockf.R --frequencyFilesList /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/frequency_files_list/n50_${RATE}.1240K.txt --popLeft1 <(cut -f 1 /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/data/n50_${RATE}/pop${POP_IDX}.group${GRP_IDX}.ind) --popRight1 pop0.x_group${GRP_IDX},pop1.x_group${GRP_IDX},pop2.x_group${GRP_IDX},pop3.x_group${GRP_IDX},pop4.x_group${GRP_IDX},pop5.x_group${GRP_IDX},pop6.x_group${GRP_IDX},pop7.x_group${GRP_IDX},pop8.x_group${GRP_IDX} --jackknife CHR --allSitesBlock --ascOutgroup Ref --outgroupMaxHetero 0 -o /mnt/archgen/users/lei_huang/RAS_tools/example_simulation/result/n50_${RATE}/f3_block.group${GRP_IDX}_pop${POP_IDX}_ind.x_group${GRP_IDX}.1240K.txt
            fi
        done
    done
done




