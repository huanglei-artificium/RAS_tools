/mnt/archgen/users/lei_huang/ancient_British/f4.Out.EUR_asc.2Test/make_freqSum.sh -> /mnt/archgen/users/lei_huang/RAS_tools/make_freqSum.sh
/mnt/archgen/users/lei_huang/ancient_British/f4.Out.EUR_asc.2Test/get.sfs2blockfsum.R -> /mnt/archgen/users/lei_huang/RAS_tools/get.sfs2blockfsum.R
/mnt/archgen/users/lei_huang/ancient_British/f4.Out.EUR_asc.2Test/get.blockfsum2blockf.R -> /mnt/archgen/users/lei_huang/RAS_tools/get.blockfsum2blockf.R


Rscript get.sfs2blockfsum.R --frequencyFilesList example/frequency_files_list.txt -n 3 --popA Outgroup --popB CEU,FIN,GBR,IBS,TSI --popBfillna Outgroup --popC example/popC.txt --jackknife CHR --ascFreqFile /mnt/archgen/public_data/HGDP+1KG_callset/stats/POS.all_filtered.AFR_mono.5EUR/count1000_floor.EUR_asc_min258.txt --outgroup AFR_all -o example/f3.AFR_all.txt