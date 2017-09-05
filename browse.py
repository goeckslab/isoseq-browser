import argparse
import json
import getGene
from bokeh.plotting import Figure
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.layouts import row, column, widgetbox
from bokeh.palettes import brewer
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, TextInput, PreText, DataTable, TableColumn, CheckboxGroup, Button, RadioButtonGroup, RadioGroup, CheckboxButtonGroup
from bokeh.models.callbacks import CustomJS
from collections import Counter

#
# Constants.
#

TITLE_FONT_SIZE = "25pt"
# color of transcripts: [reference isoorm, group1, group2...]
COLORS = brewer["Spectral"][11]
COLORS = COLORS + brewer["PuBuGn"][4]
COLORS.insert(0, '#22313F')

#
# Globals.
#

# the data used for plotting isoforms, boundaries and gene
blockDict = dict(top=[], bottom=[], left=[], right=[], exon=[],
                 start=[], end=[], chromosome=[], xs=[], ys=[], boundary=[])
tranDict = dict(top=[], bottom=[], left=[], right=[])
sourceDict = dict(xs=[], ys=[], color=[], line_alpha=[], height=[],
                  tran=[], full=[], partial=[], annot=[], start=[],
                  end=[], fileColor=[])
geneDict = dict(Gene=[], Transcripts=[])
codonDict = dict(x=[], y=[], color=[], size=[])
try:
    with open('gene.json', 'r') as f:
        data = json.load(f)
        markedDict = {'Gene': data.keys()}
        f.close()
except IOError:
    markedDict = dict()
    with open('gene.json', 'w') as f:
        json.dump(markedDict, f)
        f.close()

# update the ColumnDataSource = instant update plot
# selected exon boundaies
blockSource = ColumnDataSource(data=blockDict)
# all exon boundaries
allBlockSource = ColumnDataSource(data=blockDict)
# each transcript region
tranSource = ColumnDataSource(data=tranDict)
# the exons
source = ColumnDataSource(data=sourceDict)
# table of genes and # of clusters
geneSource = ColumnDataSource(data=geneDict)
codonSource = ColumnDataSource(data=codonDict)
markedSource = ColumnDataSource(data=markedDict)

# Create fake data source for Height and Width sliders.
slider_fake_source = ColumnDataSource(data=dict(value=[]))


# Column that holds plot.
plotColumn = column()

#
# Functions.
#


def createPlot(height=600, width=1200):
    """
    Create and return a plot for visualizing transcripts.
    """
    TOOLS = "pan, wheel_zoom, save, reset, tap"
    p = Figure(title="", y_range=[], webgl=True,
               tools=TOOLS, toolbar_location="above",
               plot_height=height, plot_width=width)
    # This causes title to overlap plot substantially:
    #p.title.text_font_size = TITLE_FONT_SIZE
    p.xgrid.grid_line_color = None               # get rid of the grid in bokeh
    p.ygrid.grid_line_color = None
    # the block of exons, there's mouse hover effect on that
    quad = p.quad(top="top", bottom="bottom", left="left", right="right",
                  source=blockSource, fill_alpha=0,
                  line_dash="dotted", line_alpha=0.4, line_color='black',
                  hover_fill_color="red", hover_alpha=0.3,
                  hover_line_color="white",
                  nonselection_fill_alpha=0, nonselection_line_alpha=0.4,
                  nonselection_line_color='black')
    # the block of each vertical transcript, each one can be selected
    p.quad(top="top", bottom="bottom", right="right", left="left",
           source=tranSource, fill_alpha=0, line_alpha=0,
           nonselection_fill_alpha=0, nonselection_line_alpha=0)
    # what exons really is
    # Cannot use line_width="height" because it is broken.
    p.multi_line(xs="xs", ys="ys", line_width=opt.height, color="color",
                 line_alpha="line_alpha", source=source)
    # the start/stop codon
    p.inverted_triangle(x="x", y="y", color="color", source=codonSource,
                        size='size', alpha=0.5)
    # mouse hover on the block
    p.add_tools(HoverTool(tooltips=[("chromosome", "@chromosome"), ("exon", "@exon"),
                ("start", "@start"), ("end", "@end")], renderers=[quad]))
    return p


def updateGene(use_saved_settings=False):
    """
    The "main" function of this app. At startup it reads the annotation and pickle file.
    When genes are changed, a new plot is created and drawn.
    """

    if opt.gene == Gene.value.strip().upper():
        geneUpdated = False
    else:
        opt.gene = Gene.value.strip().upper()                 # get the gene name from UI, pass to a global variable opt
        geneUpdated = True

    with open('gene.json', 'r') as f:
        data = json.load(f)
        if use_saved_settings and opt.gene in data.keys():
            Mark.active = [0]
            myDict = data[opt.gene]
            opt.height = myDict['height']
            Height.value = opt.height
            opt.width = myDict['width']
            Width.value = opt.width
            opt.full = myDict['full']
            Full.value = opt.full
            opt.partial = myDict['partial']
            Partial.value = opt.partial
            opt.group = myDict['group']
            Group.active = opt.group
            opt.cluster = myDict['cluster']
            Cluster.value = opt.cluster
        else:
            Mark.active = []
            opt.height = Height.value
            opt.width = Width.value
            opt.full = Full.value
            opt.partial = Partial.value
            opt.group = Group.active
            opt.cluster = Cluster.value
        f.close()

    # Clear the current plot.
    plotColumn.children = []

    # Reset the plot to blank when initial updating genes
    blockSource.data = dict(top=[], bottom=[], left=[], right=[], exon=[],
                            start=[], end=[], chromosome=[], xs=[], ys=[],
                            boundary=[])
    allBlockSource.data = dict(top=[], bottom=[], left=[], right=[], exon=[],
                               start=[], end=[], chromosome=[], xs=[], ys=[],
                               boundary=[])
    source.data = dict(xs=[], ys=[], color=[], line_alpha=[], height=[],
                       tran=[], full=[], partial=[], annot=[], start=[],
                       end=[], fileColor=[])
    codonSource.data = dict(x=[], y=[], color=[], size=[])
    matchList = Matches.value.strip().replace(' ', '').split(',')       # get the list of pickle files from UI
    opt.matches = matchList
    # load the matched isoforms from pickle file
    Console.text = 'Console:\nReading pickle file...'
    if opt.clusterDict is None:                                                 # if it's the first time to load up pickle file
        try:
            clusterDict = getGene.getMatchedIsoforms(getParams(None, matchList, None))
            opt.clusterDict = clusterDict                                       # hold pickle file dictionary in RAM
            howManyIsoforms(clusterDict, matchList)                             # find out how many isoforms for each gene
            isMatch = True                                                      # the pickle file works well
        except IOError:                                                         # if the file is not found in directory
            Console.text = 'Console:\none of the matched file \n%s is not found' % matchList
            isMatch = False
    else:
        if set(opt.clusterDict.keys()) != set(matchList):            # if the pickle files are updated, do the previous thing
            try:
                clusterDict = getGene.getMatchedIsoforms(getParams(None, matchList, None))
                opt.clusterDict = clusterDict
                howManyIsoforms(clusterDict, matchList)
                isMatch = True
            except IOError:
                Console.text = 'Console:\none of the matched file \n%s is not found' % matchList
                isMatch = False
        else:                                                       # if pickle file is not updated, do nothing
            isMatch = True
    Console.text = 'Console:\nReading annotation file...'
    if opt.annotations is None:                             # if it's the first time to load up pickle file
        try:
            opt.gtf = GTF.value.strip()                     # get the gene name from UI, pass to a global variable opt
            Annotations = getGene.getAnnotations(opt)       # get a dictionary of all transcripts in annot file, hold it in RAM
            opt.annotations = Annotations
            isAnnot = True                                  # the annotation file works well
        except IOError:
            Console.text = 'Console:\nannotations file \n%s is not found' % opt.gtf
            isAnnot = False
    else:                                                   # if the annotation files are updated, do the previous thing
        if opt.gtf != GTF.value.strip():
            try:
                opt.gtf = GTF.value.strip()
                Annotations = getGene.getAnnotations(opt)
                opt.annotations = Annotations
                isAnnot = True
            except IOError:
                Console.text = 'Console:\nannotations file \n%s is not found' % opt.gtf
                isAnnot = False
        else:                                               # if the pickle files are updated, do the previous thing
            isAnnot = True

    global tranNum, colorDF, chromosome, strand, opt
    tranList, exonList = selectGene(isAnnot, isMatch)                       # select transcripts by gene
    chromosome = getChromosome(tranList)                                         # find out which chromosome does the gene locate

    strand = exonList[0].strand                            # which strand does the gene locate on
    if strand == '+':                                      # if it's forward strand
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = getGene.assignBlocks(opt, exonList)       # assign each exon to a block
    else:                                                  # if it's trailing strand
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = getGene.assignBlocksReverse(opt, exonList)       # assign each exon to a block -- backwards

    getGene.findRegions(tranList)                       # determine regions occupied by each transcript
    tranNames = getGene.orderTranscripts(tranList)      # get the names of transcripts, placed them in the right order
    tranNames = getGene.reduceNameLength(tranNames)     # if the length of name is too long, reduce it

    tranNum = len(tranNames)                             # how many transcripts are there
    Console.text = 'Console:\nCreating plot...'

    # Create the plot to visualize gene's transcripts.
    height = int(opt.height) * 2 * (tranNum + 4)        # set plot height using transcript height
    width = int(opt.width)

    plot = createPlot(height=height, width=width)
    plotColumn.children= [plot]
    plot.title.text = "%s isoforms" % opt.gene         # update the title of plot

    # p.height = Height.value * 2 * (tranNum + 4)       # set the height of plot according to the length of transcripts
    plot.y_range.factors = tranNames[::-1]             # set the y axis tick to the transcripts names

    Console.text = 'Console:\nGrouping...'
    if 1 in opt.group and isMatch is True:
        if geneUpdated:
            colorDF = getGene.groupTran(tranList, exonList, 15)          # group the transcripts by similarities
    else:
        colorDF = None
    sourceDict = getExonData(exonList, colorDF)         # get the data of each isoform that can be directly used to plot
    codonDict = plotStartStop(tranList, blocks)         # get the location of start, stop codons
    codonSource.data = codonDict
    source.data = sourceDict

    # update the data used for plotting boundaries and hover block
    blockDict, tranDict = getBoundaryData(blocks, chromosome)              # get the data of each block that can be directly used to plot
    blockSource.data = blockDict
    allBlockSource.data = blockDict
    tranSource.data = tranDict
    if isAnnot is False:
        Console.text = 'Console:\nSuccess! Annotation\n file is missing.'
    elif isMatch is False:
        Console.text = 'Console:\nSuccess! Match file\n is missing.'
    else:
        Console.text = 'Console:\nSuccess!'


def updateGroup(attrname, old_num_clusters, new_num_clusters):
    """
    Update according to the change of grouping/clustering.
    """
    opt.cluster = Cluster.value
    opt.group = Group.active
    sourceDict = source.data
    if 0 in opt.group:                 # if it is told to group by files
        sourceDict['color'] = sourceDict['fileColor']
    else:
        if 1 in opt.group:            # if it is told to group by clustering
            colors = list()
            colors = [getColorFromDF(tran, colorDF, opt.cluster) for tran in sourceDict['tran']]
        else:
            colors = list()
            for i in sourceDict['annot']:
                if i is True:           # if it's annotation
                    colors.append(COLORS[0])
                else:
                    colors.append(COLORS[1])
        sourceDict['color'] = colors
    source.data = sourceDict


def updateHeightWidth(attrname, old, new):
    """
    Update plot height and width.
    NOTE: this method is not used right now because plot is created for each
    height/width change in updateGene(). When Bokeh bug is fixed, plot should
    be created once and updated using this function.
    """
    opt.height = Height.value
    opt.width = Width.value
    sourceDict = source.data
    sourceDict['height'] = [int(opt.height) for x in range(len(sourceDict['xs']))]
    source.data = sourceDict
    codonDict = codonSource.data
    codonDict['size'] = [int(opt.height) * 1.2 for x in range(len(codonDict['x']))]    # adjust the codon size accordingly
    codonSource.data = codonDict

    # Setting plot height, width is broken in Bokeh version 0.12.0, so this will not work:
    # plot.height = int(Height.value) * 2 * (tranNum + 4)        # update the height of plot according to the height of transcript in UI
    # plot.width = int(Width.value)                # update plot width according to width


def updateGeneTable(attrname, old, new):
    actived = Sort.active
    geneDict = geneSource.data
    df = pd.DataFrame()
    df['Gene'] = geneDict['Gene']
    df['Transcripts'] = geneDict['Transcripts']
    if actived == 0:
        df = df.sort_values(by='Gene', ascending=True)
    elif actived == 1:
        df = df.sort_values(by='Transcripts', ascending=False)
    geneDict['Gene'] = list(df['Gene'])
    geneDict['Transcripts'] = list(df['Transcripts'])
    geneSource.data = geneDict


# show/hide transcripts according to UI selection, implemented by changing the alpha values
def selectTran(attr, old, new):
    opt.full = Full.value
    opt.partial = Partial.value
    if tranSource.selected['1d']['indices'] == []:          # if no transcript is selected
        sourceDict = source.data
        # change alpha value accordingly
        sourceDict['line_alpha'] = [getAlpha(None, x) for x in zip(sourceDict['annot'],
                                    sourceDict['full'], sourceDict['partial'])]
        source.data = sourceDict
        blockSource.data = allBlockSource.data              # reset blocks to initial state
    else:
        index = tranSource.selected['1d']['indices'][0]     # which transcript is selected
        sourceDict = source.data
        # make unselected transcripts more transparent
        sourceDict['line_alpha'] = [getAlpha(index, x) for x in zip(sourceDict['ys'],
                                    sourceDict['annot'], sourceDict['full'],
                                    sourceDict['partial'])]
        source.data = sourceDict
        blocks = list()
        # in selected transcripts, find out the start and end position of each exons,
        # create new block object
        for i, yy in enumerate(sourceDict['ys']):
            if yy[0] == index + 1:
                if strand == '+':
                    start = sourceDict['start'][i]
                    end = sourceDict['end'][i]
                else:
                    start = sourceDict['end'][i]
                    end = sourceDict['start'][i]
                boundary = sourceDict['xs'][i][1]
                bound = getGene.Block(start, end, boundary)
                blocks.append(bound)
        bd, tr = getBoundaryData(blocks, chromosome)
        blockSource.data = bd


# change the alpha of each exon value according UI: full/partial widgets, select transcripts
# it's not very intuitive cause x is different in each if/else statement
def getAlpha(index, x):
    if index is None:           # nothing is selected
        if x[0] is False:       # not an annotation exon
            if x[1] < opt.full or x[2] < opt.partial:   # low full/partial reads support
                return 0
            else:
                return 1
        else:
            return 1
    else:                   # a transcript is selected
        if x[1] is False:   # if it's not an annotation exon
            if x[2] < opt.full or x[3] < opt.partial:       # low full/partial reads support
                return 0
            else:
                if index + 1 == x[0][0]:    # if the transcript is selected
                    return 1
                else:                       # if the transcript is not selected
                    return 0.3
        else:               # if it is a annotation exon
                if index + 1 == x[0][0]:    # the transcript is selected
                    return 1
                else:
                    return 0.3


# Select isoforms of a particular gene
def selectGene(isAnnot, isMatch):
    tranList = list()                              # list of Transcript objects
    exonList = list()                              # list of Exon objects
    if isAnnot:                                    # read the reference file
        try:
            getGene.getGeneFromAnnotation(opt, tranList, exonList)
        except RuntimeError:
            Console.text = 'Console:\n%s not found in annotation \nfile' % opt.gene
    if isMatch:                                    # read the pickle file
        try:
            getGene.getGeneFromMatches(opt, tranList, exonList)
        except RuntimeError:
            Console.text = 'Console:\n%s not found in pickle \nfile' % opt.gene
    return tranList, exonList


# get the data for plotting exons (start, end position for example)
def getExonData(exonList, colorDF):
    sourceDict = dict(xs=[], ys=[], color=[], line_alpha=[], height=[],
                      tran=[], full=[], partial=[], annot=[], start=[],
                      end=[], fileColor=[])
    columns = ['xs', 'ys', 'color', 'start', 'end', 'tran', 'full',
               'partial', 'annot', 'fileColor']
    num_clusters = opt.cluster
    for myExon in exonList:
        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart

        if 0 in opt.group:       # if group by files, pass
            color = COLORS[myExon.tran.source[0]]
        else:
            if colorDF is not None:
                # Get transcript color based on grouping.
                color = getColorFromDF(myExon.tran.name, colorDF, num_clusters)
            else:          # if the grouping effect is off, paint default color
                if myExon.tran.annot:
                    color = COLORS[0]
                else:
                    color = COLORS[1]
        xs = (adjStart, adjStart + exonSize)
        ys = (tranNum - (myExon.tran.tranIx), tranNum - (myExon.tran.tranIx))
        values = [xs, ys, color, myExon.start, myExon.end,
                  myExon.tran.name, myExon.tran.full, myExon.tran.partial,
                  myExon.tran.annot, COLORS[myExon.tran.source[0]]]

        for ix, col in enumerate(columns):
            sourceDict[columns[ix]].append(values[ix])

    sourceDict['line_alpha'] = [1 for x in range(len(sourceDict['xs']))]
    sourceDict['height'] = [int(opt.height) for x in range(len(sourceDict['xs']))]
    return sourceDict


def getColorFromDF(transcript_name, colorDF, num_clusters):
    """
    Get transcript color based on number of clusters.
    """
    if transcript_name not in list(colorDF.name):
        color = COLORS[0]
    else:
        row = colorDF.loc[colorDF['name'] == transcript_name]    # find out which transcript it is, and what group it belongs
        groupName = 'group%d' % num_clusters          # how many groups are there
        try:
            group = row[groupName].values[0]
            color = COLORS[group + 1]
        except (ValueError, KeyError):          # if the input groups are more than total number of transcripts
            color = COLORS[1]
    return color


# find out the position of boundaries
def getBoundaryData(blocks, chromosome):
    blockDict = dict(top=[], bottom=[], left=[], right=[], exon=[],
                     start=[], end=[], chromosome=[], xs=[], ys=[],
                     boundary=[])
    tranDict = dict(top=[], bottom=[], left=[], right=[])
    blockDict['xs'] = [(0, 0)]
    blockDict['ys'] = [(0, tranNum + 1)]
    columns = ['boundary', 'left', 'right', 'xs', 'start', 'end', 'exon']
    exonCounter = 1
    numberOfBlocks = len(blocks)
    for bound in blocks:      # infomation for the mouse hover effect on blocks
        if strand == '+':
            values = [bound.boundary, bound.boundary + bound.start - bound.end,
                      bound.boundary, (bound.boundary, bound.boundary),
                      bound.start, bound.end, exonCounter]
        else:
            values = [bound.boundary, bound.boundary - bound.start + bound.end,
                      bound.boundary, (bound.boundary, bound.boundary),
                      bound.start, bound.end, exonCounter]
        blockDict['top'] = [(tranNum + 1) for x in range(numberOfBlocks)]
        blockDict['bottom'] = [0 for x in range(numberOfBlocks)]
        blockDict['chromosome'] = [chromosome for x in range(numberOfBlocks)]
        exonCounter += 1
        counter = 0
        while counter < len(columns):
            blockDict[columns[counter]].append(values[counter])
            counter += 1
    right = blockDict['right']
    blockDict['ys'] = [(0, tranNum + 1) for x in range(numberOfBlocks + 1)]

    # put the region of each transcript into a block
    tranDict['top'] = [x + 1.5 for x in range(tranNum)]
    tranDict['bottom'] = [x + 0.5 for x in range(tranNum)]
    tranDict['left'] = [0 for x in range(tranNum)]
    tranDict['right'] = [max(right) for x in range(tranNum)]
    return blockDict, tranDict


# find out the chromosome that isosoforms locate on, find by matched isoform
def getChromosome(tranList):
    chromosome = None
    for tran in tranList:
        if tran.annot is False:               # find it in the matched isoforms
            chromosome = tran.chr
            break
    return chromosome


def howManyIsoforms(clusterDict, matchList):
    allGenes = Counter()                                        # create a counter hastable(dictionary) object
    for matchFile in matchList:
        geneDict = dict()
        myDict = clusterDict[matchFile].getGeneDict()
        for key, val in myDict.iteritems():
            geneDict.setdefault(key, len(val))                          # how many isoforms for each gene
        geneDict = Counter(geneDict)
        allGenes = allGenes + geneDict                                  # combine every match files
    df = pd.DataFrame()
    df['Gene'] = allGenes.keys()
    df['Transcripts'] = allGenes.values()
    df = df.sort_values(by='Transcripts', ascending=False)
    geneSource.data = dict(Gene=list(df['Gene']), Transcripts=list(df['Transcripts']))


# save the transcripts to .fasta file, the function is copied from MatchAnnot
def saveFasta(attrname, old, new):
    Console.text = 'Console:\nSaving...'
    opt.fasta = Save.value.strip()
    tranList = list()
    exonList = list()
    getGene.getGeneFromMatches(opt, tranList, exonList)
    opt.fasta = None
    Console.text = 'Console:\nSuccessfully saved'


def markGene(attrname, old, new):
    if 0 in Mark.active:
        with open('gene.json', 'r') as f:
            data = json.load(f)
            f.close()
        with open('gene.json', 'w+') as f:
            paramDict = {'height': opt.height, 'width': opt.width, 'full': opt.full,
                         'partial': opt.partial, 'cluster': opt.cluster, 'group': opt.group}
            if opt.gene not in data.keys():
                data[opt.gene] = paramDict
            f.write(json.dumps(data))
            markedSource.data = dict(Gene=data.keys())
            f.truncate()
            f.close()
    else:
        with open('gene.json', 'r') as f:
            data = json.load(f)
            f.close()
        with open('gene.json', 'w+') as f:
            if opt.gene in data:
                data.pop(opt.gene, None)
            f.write(json.dumps(data))
            markedSource.data = dict(Gene=data.keys())
            f.truncate()
            f.close()


def plotStartStop(tranList, blocks):
    '''Add start/stop codons to plot.'''
    codonDict = dict(x=[], y=[], color=[])
    for tran in tranList:
        if tran.annot:                             # only annotations know about start/stops
            if hasattr(tran, 'startcodon'):
                codonDict['color'].append('green')
                xPos = findCodon(tran.startcodon, blocks)
                codonDict['x'].append(xPos)
                codonDict['y'].append(tranNum - tran.tranIx)
            if hasattr(tran, 'stopcodon'):
                codonDict['color'].append('red')
                xPos = findCodon(tran.stopcodon, blocks)
                codonDict['x'].append(xPos)
                codonDict['y'].append(tranNum - tran.tranIx)
    codonDict['size'] = [int(opt.height) * 1.2 for x in range(len(codonDict['x']))]
    return codonDict


def findCodon(posit, blocks):
    '''Add a codon mark to the plot.'''
    for blk in blocks:
        if blk.start <= posit and blk.end >= posit or \
                blk.start >= posit and blk.end <= posit:      # check in both strand directions
            xPos = blk.boundary - abs(blk.end - posit)
    return xPos

#
# Classes.
#


# a class containing all the input parameters
class getParams(object):
    def __init__(self, gtf, matches, gene, format="standard", fasta=None,
                 annotations=None, clusterDict=None, height=None, width=None,
                 full=None, partial=None, group=None, cluster=None):
        self.gtf = gtf                              # reference genome file
        self.matches = matches                      # list of matched files
        self.gene = gene                            # which gene to load
        self.format = format                        # format is a keyword
        self.fasta = fasta                          # output isoforms to fasta
        self.annotations = annotations              # hold annotation dictionary in RAM
        self.clusterDict = clusterDict              # hold isoforms information in RAM
        self.height = height
        self.width = width
        self.full = full
        self.partial = partial
        self.group = group
        self.cluster = cluster


#
# Application code.
#

# Read command line arguments and use to set widget defaults.
parser = argparse.ArgumentParser(description='Visual analytics for PacBio data.')
parser.add_argument('--input', dest='input_file', help='Input file (pickle)')
parser.add_argument('--anno', dest='anno_file', help='Annotation file (gtf)')
args, unknown = parser.parse_known_args()
input_file = args.input_file or "matches.pickle"
anno_file = args.anno_file or "gencode.vM9.annotation.gtf"

# Create all widgets.
GTF = TextInput(title="Annotation file", value=anno_file)
Format = TextInput(title="Annotation file format, standard is gtf", value="standard")
Matches = TextInput(title="MatchAnnot pickle files (ex: a.pickle, b.pickle)", value=input_file)
Gene = TextInput(title="Gene to visualize", value="BRCA1")
Full = Slider(title="Full reads support threshold",
              value=0, start=0, end=30, step=1.0)
Partial = Slider(title="Partial reads support threshold",
                 value=0, start=0, end=50, step=1.0)
Group = CheckboxGroup(labels=["Group by file", "Group by similarity"],
                      active=[1])
Cluster = Slider(title="Number of isoform groups",
                 value=3, start=1, end=15, step=1.0)
Height = Slider(title="Transcript height", value=10, start=5, end=30, step=1)
Width = Slider(title="Plot width", value=600, start=400, end=1500, step=50)
Save = TextInput(title="Enter the folder name to save data in Fasta", value=None)
button = Button(label='GO', button_type="success")
Sort = RadioButtonGroup(labels=["Rank by Gene", "Rank by Transcripts"], active=1)
Mark = CheckboxButtonGroup(labels=["Save gene"], active=[])

opt = getParams(None, [], None, format=None)    # a object that contains all the inputs options for read data

# the console box
Console = PreText(text='Console:\nStart visualize by entering \nannotations, pickle file and\n gene. Press Enter to submit.\n', height=70)


# a table of with all the genes in the match files, and how many isoforms in each gene
geneColumns = [TableColumn(field="Gene", title="Gene"),
               TableColumn(field="Transcripts", title="Isoforms")]
geneCountTable = DataTable(source=geneSource, columns=geneColumns, sortable=False,
                           row_headers=False, width=280)
markedColumns = [TableColumn(field="Gene", title="Saved genes")]
markedGeneTable = DataTable(source=markedSource, columns=markedColumns, sortable=False, width=280, row_headers=False)

# make changes to the plot when widgets are updated
button.on_click(updateGene)
Mark.on_change('active', markGene)
Full.on_change('value', selectTran)
Partial.on_change('value', selectTran)
Cluster.on_change('value', updateGroup)
Group.on_change('active', updateGroup)
Save.on_change('value', saveFasta)
tranSource.on_change('selected', selectTran)
Sort.on_change('active', updateGeneTable)

# Add mouseup callback on sliders.
slider_fake_source.on_change('data', lambda attr, old, new: updateGene())
for slider in [Height, Width]:
    slider.callback_policy = "mouseup"
    slider.callback = CustomJS(args=dict(source=slider_fake_source), code="""
        source.data = { value: [cb_obj.value] }
    """)

# Add handlers for selecting genes from tables. Handlers update the Gene textinput
# and updates the plot.
def add_selected_handler(table, use_saved_settings):
    data_source = table.source

    # Callback updates Gene value with selected gene.
    def internal_callback(attr, old, new):
        selected_index = new["1d"]["indices"][0]
        Gene.value = data_source.data['Gene'][selected_index]
        updateGene(use_saved_settings)
    data_source.on_change("selected", internal_callback)

add_selected_handler(geneCountTable, False)
add_selected_handler(markedGeneTable, True)


# Layout interface.
inputs_and_outputs = [Console, GTF, Matches, Format, Save]
plot_controls = [Gene, button, Group, Cluster, Full, Partial, Height, Width, Sort, geneCountTable, Mark, markedGeneTable]

curdoc().add_root(row( row(inputs_and_outputs), row(widgetbox(plot_controls), plotColumn) ) )

curdoc().add_root(slider_fake_source)
curdoc().title = "Iso-Seq Browser"
