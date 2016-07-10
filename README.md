# Isoseq-browser

![alt text](images/BRCA1.png)

Isoseq-browser is a interactive visualization tool for RNA isoforms. It is build on [bokeh](https://github.com/bokeh/bokeh) visualization library of Python. Isoseq-browser is designed for Pacific Bioscience [Isoseq](http://www.pacb.com/blog/intro-to-iso-seq-method-full-leng/). Isoseq-browser is build on [MatchAnnot](https://github.com/TomSkelly/MatchAnnot), which is match isoforms to annotation genome.

# Install
* Clone the repositories or download zip
* Make environment:

```
    make env Makefile
```

* Download the example data (optional):

```
    make download Makefile
```

* Run Isoseq-browser:

```
    make run
```

# How to run Isoseq-browser
* Run MatchAnnot, the running instructions are [here](https://github.com/TomSkelly/MatchAnnot/wiki/How-to-Run-matchAnnot). Isoseq-browser requires the pickle file of MatchAnnot.
* Run Isoseq-browser:

```
    make run
```

* Enter the annotation file and matched pickle file
* Enter the gene to start visualize
* Use the other widgets for more options:

# Tips
* It will take longer time to load the first gene, but it take shorter time when change genes.
* For the grouping isoforms. It will take a while to group more than 40 clusters. The grouping algorithm is originally come from [Isoview](https://github.com/JMF47/IsoView).
