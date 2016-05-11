# For setting up a conda environment.
ENV_NAME=ib_env
ACTIVATE_ENV=bash -c "source activate $(ENV_NAME)"

run: env gencode.vM9.annotation.gtf
	# python matchAnnot.py --gtf gencode.vM9.annotation.gtf --outpickle True example.sam > example.pickle
	PYTHONPATH=./dep bokeh serve --show browse.py

# Set up the environment and dependencies.
env:
	mkdir dep && git clone https://github.com/TomSkelly/MatchAnnot dep
	conda create -n $(ENV_NAME) -y python && conda install -y -c bokeh scikit-learn
	touch env

gencode.vM9.annotation.gtf:
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gtf.gz
	gunzip gencode.vM9.annotation.gtf.gz


clean:
	rm -rf dep
	rm -f gencode.vM9.annotation.gtf.gz
	conda remove --name $(ENV_NAME) --all -y
