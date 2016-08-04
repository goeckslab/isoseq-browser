SHELL=bash

#
# Parameters.
#

MATCHES_PICKLE_FROM_SRC=mcf7_matchAnnot_results_src.pickle
MATCHES_PICKLE_FROM_DOWNLOAD=mcf7_matchAnnot_results_download.pickle
MATCHES_INPUT=$(MATCHES_PICKLE_FROM_DOWNLOAD)
ANNOTATION_VERSION=25
ANNOTATION_GTF=gencode.v$(ANNOTATION_VERSION).annotation.gtf

# For setting up a conda environment.
ENV_NAME=ib_env
ACTIVATE_ENV=source activate $(ENV_NAME)

run: env $(MATCHES_INPUT) $(ANNOTATION_GTF)
	$(ACTIVATE_ENV) && PYTHONPATH=./dep:. bokeh serve --show browse.py --args --input $(MATCHES_INPUT) --anno $(ANNOTATION_GTF)

# Download and unzip GENCODE annotation.
$(ANNOTATION_GTF):
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_$(ANNOTATION_VERSION)/$(ANNOTATION_GTF).gz
	gunzip $(ANNOTATION_GTF).gz

# Set up the environment and dependencies.
env:
	mkdir dep && git clone https://github.com/TomSkelly/MatchAnnot dep
	conda create -n $(ENV_NAME) -y python
	$(ACTIVATE_ENV) && conda install -y pandas bokeh=0.12.0 scikit-learn
	touch env

# Download precomputed MatchAnnot pickle file.
$(MATCHES_PICKLE_FROM_DOWNLOAD):
	wget -O $(MATCHES_PICKLE_FROM_DOWNLOAD) http://goeckslab.org/files/mcf7_matchAnnot_results.pickle

# Remove and clean up everything.
clean:
	rm -f mcf7_matchAnnot_results.*
	rm -f env
	rm -rf dep
	conda remove --name $(ENV_NAME) --all -y

# Create MatchAnnot pickle file needed by browser from original sources.
# NOTE: this is computationally intensive.
$(MATCHES_PICKLE_FROM_SRC): $(ANNOTATION_GTF)
	# Download PacBio MCF7 dataset.
	wget http://datasets.pacb.com.s3.amazonaws.com/2015/IsoSeqHumanMCF7Transcriptome/IsoSeq_MCF7_2015edition_polished.unimapped.fasta

	# Install STAR and SAMtools
	$(ACTIVATE_ENV) && conda install -c bioconda -y star samtools

	# Download genome.
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_$(ANNOTATION_VERSION)/GRCh38.p7.genome.fa.gz
	gunzip GRCh38.p7.genome.fa.gz

	# Create STAR genome database.
	mkdir -p gencode/human/GRCh38.p7/star/
	$(ACTIVATE_ENV) && STAR --runThreadN 15 --runMode genomeGenerate \
		--genomeDir gencode/human/GRCh38.p7/star \
		--genomeFastaFiles GRCh38.p7.genome.fa \
		--sjdbGTFfile $(ANNOTATION_GTF) \
		--sjdbOverhang 100

	# Align reads using STAR long:
	$(ACTIVATE_ENV) && STARlong --runThreadN 15 \
			--genomeDir gencode/human/GRCh38.p7/star \
			--runMode alignReads \
			--outSAMattributes NH HI NM MD \
			--readNameSeparator space \
			--outFilterMultimapScoreRange 1 \
			--outFilterMismatchNmax 2000 \
			--scoreGapNoncan -20 \
			--scoreGapGCAG -4 \
			--scoreGapATAC -8 \
			--scoreDelOpen -1 \
			--scoreDelBase -1 \
			--scoreInsOpen -1 \
			--scoreInsBase -1 \
			--alignEndsType Local \
			--seedSearchStartLmax 50 \
			--seedPerReadNmax 100000 \
			--seedPerWindowNmax 1000 \
			--alignTranscriptsPerReadNmax 100000 \
			--alignTranscriptsPerWindowNmax 10000 \
			--readFilesIn IsoSeq_MCF7_2015edition_polished.unimapped.fasta

	# Prepartion for MatchAnnot: (1) sort aligned reads and (2) get list of valid
	# contigs, which are contigs that have annotations. MatchAnnot produces an error
	# for reads aligned to contigs that do not have annotations.
	$(ACTIVATE_ENV) && samtools view -b Aligned.out.sam | samtools sort - > mcf7_aligned.bam
	grep -v ^# $(ANNOTATION_GTF) | cut -f1 | uniq > valid_matchannot_contigs.txt

	# Use MatchAnnot to combine aligned reads and annotation.
	$(ACTIVATE_ENV) && samtools view mcf7_aligned.bam | \
	python filter_sam_by_contigname.py valid_matchannot_contigs.txt | \
	dep/matchAnnot.py --gtf $(ANNOTATION_GTF) --outpickle mcf7_matchAnnot_results_src.pickle > mcf7_matchAnnot_results_src.txt
