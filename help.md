
# IsoSeq Browser Help.

|   |   |   |   |   |
|---|---|---|---|---|
|   |   |   |   |   |
|   |   |   |   |   |
|   |   |   |   |   |

* `annotation` : annotations file, in format specified by --format, Reload page to update
* `format`:
* `matches`:

paramDict = dict(Parameter=['annotation', 'format', 'matches', 'gene', 'full', 'alpha',
                            'partial', 'group', '# of groups', 'fasta'],
                 Description=['',
                              'format of annotation file: standard, alt, pickle (default=standard: gtf)',
                              'pickle file from matchAnnot.py, if there are multiple files, seperate them with comma no space, Reload page to update',
                              'gene to plot (required)',
                              'add alpha channel to transcripts without low full/partial support, change full/partial to update',
                              'full support threshold, work together with alpha',
                              'partial support threshold, work together with alpha',
                              'group the trascript on similarity',
                              'assign transcripts into how many groups',
                              'out put the .fasta files to a folder, enter the folder name'],)
