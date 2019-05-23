featureCounts -t exon -D 1500 -p --primary -T 3 -s 1 -a mm9_with_ERCC.gff -o Aging_Heart_counts_genes.txt \
	3m4_heart_ATCACGAligned.out.bam 3m5_heart_CGATGTAligned.out.bam 3m6_heart_TTAGGCAligned.out.bam \
	12m4_heart_TGACCAAligned.out.bam 12m5_heart_ACAGTGAligned.out.bam 12m6_heart_GCCAATAligned.out.bam \
	29m4_heart_CAGATCAligned.out.bam 29m5_heart_ACTTGAAligned.out.bam 29m6_heart_GATCAGAligned.out.bam
