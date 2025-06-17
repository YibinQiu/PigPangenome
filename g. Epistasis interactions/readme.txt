01.adjust_phe.sbatch    #adjust phenotype (top five principal components, slaughter batch, sex, pre-slaughter weight were used)
02.epistasis.sbatch          #Epistasis interactions were tested using PLINK v1.9 between GWAS hit (p-value < 0.001)
03.ld_epi.sbatch     #Caculated ld (window = 2 MB) 
04.parse_ld_epi_exclude.sbatch     #parse each epi pair with ld
05.pick_tag_epi.sbatch            #defined interaction block, with the most significant pair chosen as the representative
06.adjust_tag_epi.sbatch             #Bonferroni correction for representative pairs
adjust_phe.R            #adjust phenotype (a Rscript used to 01.adjust_phe.sbatch)
adjust_tag_epi.R             #Bonferroni correction for representative pairs (a Rscript used to 06.adjust_tag_epi.sbatch)
parse_ld_epi_exclude.R             #parse each epi pair with ld (a Rscript used to 04.parse_ld_epi_exclude.sbatch)
pick_tag_epi.py             #chosen representative pair (a python script used to 05.pick_tag_epi.sbatch)