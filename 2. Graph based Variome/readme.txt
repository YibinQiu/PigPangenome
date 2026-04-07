01.PanSV_genotyping.sbatch    #vg giraffe SV genotyping
02.Multiallele2biallele_and_SelectPASS.sbatch          #Filter SVs
03.Merge_PanSVSample.sbatch     #Merge individual-based SVs to population-based SVs data
04.PanSV_QC.sbatch     #SVs quality control
05.PanSNP_PanIndel_genotyping.sbatch             #SNP and Indel calling followed by https://github.com/qgg-lab/swim-public
06.PanSNP_PanIndel_QC.sbatch             #SNP and Indel quality control
07.PigPan_PCA.sbatch             #SNP and Indel, SVs PCA
08.PigPan_NJtree.sbatch             #SNP and Indel, SVs PCA NJ tree