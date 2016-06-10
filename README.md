# Isoseq-browser
Isoseq-browser is an gene isoform visualization tool. It is build on [bokeh](https://github.com/bokeh/bokeh) visualization library of Python. Isoseq-browser is designed for Pacific Bioscience [Isoseq](http://www.pacb.com/blog/intro-to-iso-seq-method-full-leng/). Isoseq-browser is build on [MatchAnnot](https://github.com/TomSkelly/MatchAnnot), which is match isoforms to annotation genome. 

# Install
1. Clone the repositories or download zip
2. Make environment:

```
    make env Makefile
```

3. Download the example data (optional):

```
    make download Makefile
```

4. Run Isoseq-browser:

```
    make run
```

# How to run Isoseq-browser
1. Run MatchAnnot, the running instructions are [here](https://github.com/TomSkelly/MatchAnnot/wiki/How-to-Run-matchAnnot). Isoseq-browser requires the pickle file of MatchAnnot.
2. Run Isoseq-browser:

```
    make run
```

2. Enter the annotation file and matched pickle file
3. Enter the gene to start visualize
4. Use the other widgets for more options:
    
# Tips
1. It will take longer time to load the first gene, but it take shorter time to 
2. For the grouping isoforms. It will take a while to group more than 40 clusters. The grouping algorithm is originally come from [Isoview].
