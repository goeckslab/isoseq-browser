SHELL=bash

# For setting up a conda environment.
ENV_NAME=ib_env
ACTIVATE_ENV=source activate $(ENV_NAME)

run: env gencode.vM9.annotation.gtf
	# python matchAnnot.py --gtf gencode.vM9.annotation.gtf --outpickle True example.sam > example.pickle
	$(ACTIVATE_ENV) && PYTHONPATH=./dep:. bokeh serve --show browse.py

# Set up the environment and dependencies.
env:
	mkdir dep && git clone https://github.com/TomSkelly/MatchAnnot dep
	conda create -n $(ENV_NAME) -y python
	$(ACTIVATE_ENV) && conda install -y -c bokeh scikit-learn
	touch env

# Download and unzip mm10 (?) annotation.
gencode.vM9.annotation.gtf:
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gtf.gz
	gunzip gencode.vM9.annotation.gtf.gz

# Download and unzip hg19 annotation.
gencode.v24.chr_patch_hapl_scaff.annotation.gtf:
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz
	gunzip gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz

clean:
	rm -f env
	rm -rf dep
	rm -f gencode.vM9.annotation.gtf.gz
	conda remove --name $(ENV_NAME) --all -y
