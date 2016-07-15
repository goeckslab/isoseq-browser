SHELL=bash

# Input files.
MATCHES_PICKLE=matches.pickle
ANNOTATION=gencode.v24.annotation.gtf

# For setting up a conda environment.
ENV_NAME=ib_env
ACTIVATE_ENV=source activate $(ENV_NAME)

run: env $(ANNOTATION)
	# python matchAnnot.py --gtf gencode.vM9.annotation.gtf --outpickle True example.sam > example.pickle
	$(ACTIVATE_ENV) && PYTHONPATH=./dep:. bokeh serve --show browse.py --args --input $(MATCHES_PICKLE) --anno $(ANNOTATION)

# Set up the environment and dependencies.
env:
	mkdir dep && git clone https://github.com/TomSkelly/MatchAnnot dep
	conda create -n $(ENV_NAME) -y python
	$(ACTIVATE_ENV) && conda install -y pandas bokeh=0.12.0 scikit-learn
	touch env

# Download and unzip annotation from GENCODE.
$(ANNOTATION):
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/$(ANNOTATION).gz
	gunzip $(ANNOTATION).gz

clean:
	rm -f env
	rm -rf dep
	conda remove --name $(ENV_NAME) --all -y
