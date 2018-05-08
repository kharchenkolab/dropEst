all: transcripts folders

transcripts: SCG_71_C2_S2_R1_001.fastq.gz SCG_71_C2_S2_R2_001.fastq.gz SCG_71_C2_S2_R3_001.fastq.gz
genes: mm10_genes.gtf
folders:
	mkdir -p 01_dropTag 02_alignment 03_dropEst

SCG_71_C2_S2_R1_001.fastq.gz SCG_71_C2_S2_R2_001.fastq.gz SCG_71_C2_S2_R3_001.fastq.gz:
	wget http://pklab.med.harvard.edu/viktor/dropest_paper/scg71_demo/SCG_71.tar.gz
	tar xf SCG_71.tar.gz
	rm SCG_71.tar.gz

mm10_genes.gtf:
	wget http://pklab.med.harvard.edu/viktor/dropest_paper/scg71_demo/mm10_genes.gtf.gz
	gunzip mm10_genes.gtf.gz

results:
	wget http://pklab.med.harvard.edu/viktor/dropest_paper/scg71_demo/results.tar.gz
	mkdir results
	tar xf results.tar.gz -C results
	rm results.tar.gz

