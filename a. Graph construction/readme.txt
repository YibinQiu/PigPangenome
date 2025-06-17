01.ONT_mapping.sbatch    #ONT reads mapped to the W64 Yorkshire reference
02.ONT_SVdetection.sbatch          #SVs were called using Sniffles
03.genome_alignment.sbatch     #Genome alignment was performed using MuMmer
04.genome_SVdetection.sbatch     #SVs were called using SVMU
05.svmu2vcf             #SVMU format convert to vcf
06.filter_ONT.sbatch             #ONT SVs filter
07.filter_Assemble.sbatch             #Genome SVs filter
08.merge_ONT_merge_Assemble.sbatch             #Merge SVs
09.bench_ONT_vs_Chromosome_Scaffold_Contig.sbatch             #Find overlaped SVs
10.remove_Chromosome_Scaffold_Contig_overlapSV.sbatch             #Remove overlap SVs from genome origin
11.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV.sbatch             #Merge all SVs
12.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV_all_0.8.sbatch             #Trvari collapse SV
13.vg_autoindex_0.8.sbatch             #SV graph construction
svmu2vcf.R             #SVMU format convert to vcf (a Rscript used to 05.svmu2vcf)