# Iso-Seq Browser Help.

## Visualization parameters

| Parameter  |  Description  |
|---|---|
| `annotation`  | Annotations file, in format specified by --format. Reload page to update e.g.*example.gtf*  |
| `format`  | Format of annotation file: standard (gtf), alt, pickle  |
| `matches`  | Pickle file from [MatchAnnot](https://github.com/TomSkelly/MatchAnnot). For multiple files, separate them with comma. Reload page to update. e.g. *match1.pickle,match2.pickle* |
| `Isoform height` | Height of each isoform/transcript |
| `Plot width` | Width of the plot |
| `full`  | Full-length read threshold, transcripts with lower full supports will not be displayed   |
| `partial` | Partial-length read threshold, transcripts with lower partial supports will not be displayed  |
| `Group by file` | Group transcripts by different files (when there are more than one matches file) |
| `group`  | Group the transcripts by similarity (using K-Means algorithm) |
| `number of groups`  | Assign transcripts into how many groups  |
| `fasta`  | Folder name for fasta output files of exported data  |


## Gene table parameters

| Parameter  |  Description  |
|---|---|
| `Gene`  | The gene to visualize (required)  |
| `Confirm button` | update the visualization |
| `Enter from box` / `Select from geneTable` | Isoseq-browser allows two way to input gene. One way is to type the gene name into the text box, the other is to select gene from the generated geneTable, which has information of all genes and count of transcripts for that gene. |

## Bokeh plot tools located to the upper right of the visualization

| Tool  |  Description  |
|---|---|
| `pan`  | move around the plot  |
| `wheel zoom` | zoom in/out using mouse wheel or touchpad scrolling |
| `tap` | select one transcript to see its exons |
| `save` | save the plot as image |
| `reset` | reset the plot by removing changes made using Bokeh tools such as zoom and move |
| `hover tool` | hover on exon will display more information, it can be turned off |
