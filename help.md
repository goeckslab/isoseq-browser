# Iso-Seq Browser Help.

## Visualization parameters

| Parameter  |  Description  |
|---|---|
| `Annotation`  | Annotations file, in format specified by --format. Reload page to update e.g.*example.gtf*  |
| `Matches`  | Pickle file from [MatchAnnot](https://github.com/TomSkelly/MatchAnnot). For multiple files, separate them with comma. Reload page to update. e.g. *match1.pickle,match2.pickle* |
| `Format`  | Format of annotation file: standard (gtf), alt, pickle  |
| `Fasta`  | Folder name for fasta output files of exported data  |
| `Transcript height` | Height of each isoform/transcript |
| `Plot width` | Width of the plot |
| `Full`  | Full-length read threshold, transcripts with lower full supports will not be displayed   |
| `Partial` | Partial-length read threshold, transcripts with lower partial supports will not be displayed  |
| `Group by file` | Group transcripts by different files (when there are more than one matches file) |
| `Group by similarity`  | Group the transcripts by similarity (using K-Means algorithm) |
| `number of groups`  | Assign transcripts into how many groups  |


## Gene table parameters

| Parameter  |  Description  |
|---|---|
| `Enter from box` / `Select from geneTable` / `Select from marked genes` | Isoseq-browser allows three way to input gene. One way is to type the gene name into the text box, others are to select gene from the generated geneTable or marked geneTable, which has information of all genes and count of transcripts for that gene. |
| `Gene to visualize`  | The gene to visualize (required)  |
| `Go button` | update the visualization |
| `Rank Transcript` | Sort the geneTable |
| `Mark` | `Mark genes and save their input parameters for future usage` |


## Bokeh plot tools located to the upper right of the visualization

| Tool  |  Description  |
|---|---|
| `pan`  | move around the plot  |
| `wheel zoom` | zoom in/out using mouse wheel or touchpad scrolling |
| `tap` | select one transcript to see its exons |
| `save` | save the plot as image |
| `reset` | reset the plot by removing changes made using Bokeh tools such as zoom and move |
| `hover tool` | hover on exon will display more information, it can be turned off |
