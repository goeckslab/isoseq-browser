
# IsoSeq Browser Help.

| parameter  |  Description  |
|---|---|
| `annotation`  | annotations file, in format specified by --format, Reload page to update e.g.*example.gtf*  |
| `format`  | format of annotation file: standard, alt, pickle  |
| `matches`  | pickle file from matchAnnot.py, if there are multiple files, seperate them with comma no space, Reload page to update e.g. *match1.pickle,match2.pickle* |
| `Transcripts height` | height of each transcript |
| `Plot width` | width of the plot |
| `full`  | full-length read threshold, transcripts with lower full supports will not be displayed   |
| `partial` | partial-length read threshold, transcripts with lower partial supports will not be displayed  |
| `Group by file` | group transcripts by different files (when there are more than one matches file) |
| `group`  | group the transcripts by similarity (using K-Means algorithm) |
| `number of groups`  | assign transcripts into how many groups  |
| `fasta`  | out put the .fasta files to a folder, enter the folder name  |

| GeneTable parameters  |  Description  |
|---|---|
| `Select gene to visualize`  | plot which gene (required)  |
| `Confirm button` | update the plot |
| `Enter from box` / `Select from geneTable` | Isoseq-browser allows two way to input gene. One way is to type the gene name into the text box, the other is to select gene from the generated geneTable, which has information of all genes and count of transcripts for that gene. |

| Bokeh plot tools  |  Description  |
|---|---|
| `pan`  | move around the plot  |
| `wheel zoom` | zoom in/out |
| `tap` | select one transcript to see it's exons |
| `save` | save the plot as image |
| `reset` | reset the plot, get rid of bokeh tools made changes (zoom, move for instance) |
| `hover tool` | hover on exon will display more information, it can be turned off |
