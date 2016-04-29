env:
	mkdir dep
	git clone https://github.com/TomSkelly/MatchAnnot dep
	conda install -c bokeh scikit-learn

download:
	wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gtf.gz
	tar -zxvf gencode.vM9.annotation.gtf.gz

run: env
	python matchAnnot.py --gtf gencode.vM9.annotation.gtf --outpickle True example.sam > example.pickle
	bokeh serve --show browse.py

clean:

	rm gencode.vM9.annotation.gtf.gz
