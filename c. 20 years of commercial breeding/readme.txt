01.PigPan_Fst.sbatch    #Multi-generation Fst
02.Fst_filter.sbatch          #Retain Fst > 0.1
03.Caculate_Fst_Ref_freq.sbatch     #Caculate reference allele frequency (Fst > 0.1)
04.MultiGen_SimGeneFlow.sbatch     #Simulate gene flow (S21 Duroc and W51 Landrace)
05.MultiGen_SimGeneDrift.sbatch             #Simulate genetic drift (S21 Duroc and W51 Landrace)
06.SelectionPrefernece_Fst.sbatch             #Filter variants and classify into four trajectory patterns
07.QTLAnnotation_Fst.sbatch             #QTL enrichment for trajectory patterns of artificial selection variants
08.QTLDirection_Duroc_Landrace.R             #The proportion of variants with identical trajectory patterns was calculated, and defined QTL convergent and divergent direction
09.Genetic_gain_per_year.R             #Caculated genetic gain per year (S21 Duroc and W51 Landrace)
10.XPCLR_WENS_vs_USA.sbatch             #Differentiation between Chinese and American populations (XP-CLR)
11.Fst_WENS_vs_USA.sbatch             #Differentiation between Chinese and American populations (Fst)
12.QTLAnnotation_WENS_vs_USA.sbatch             #QTL enrichment analysis for selection windows
Fst_filter.R            #A Rscript for Fst > 0.1
MultiGen_SimGeneDrift_S21.prm             #Simulate genetic drift file (S21 Duroc) (a prm file used to 05.MultiGen_SimGeneDrift.sbatch)
MultiGen_SimGeneDrift_W51.prm             #Simulate genetic drift file (W51 Landrace) (a prm file used to 05.MultiGen_SimGeneDrift.sbatch)
MultiGen_SimGeneFlow_S21.prm            #Simulate genetic flow file (S21 Duroc) (a prm file used to 04.MultiGen_SimGeneFlow.sbatch)
MultiGen_SimGeneFlow_W51.prm             #Simulate genetic flow file (S21 Landrace) (a prm file used to 04.MultiGen_SimGeneFlow.sbatch)
QTLAnnotation_Fst.R             #QTL enrichment (a Rscript used to 07.QTLAnnotation_Fst.sbatch)
QTLAnnotation_WENS_vs_USA.R             #QTL enrichment (a Rscript used to 12.QTLAnnotation_WENS_vs_USA.sbatch)
SelectionPrefernece_Fst.R             #Filter variants and classify into four trajectory patterns (a Rscript used to 06.SelectionPrefernece_Fst.sbatch)