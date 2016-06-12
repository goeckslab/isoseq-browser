from getGene import *
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HoverTool, HBox, VBoxForm
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, Select, TextInput, PreText, DataTable, TableColumn
from collections import Counter

PLOT_WIDTH = 1200                                           # the width of plott
TITLE_FONT_SIZE = "25pt"                                    # font size of plot title
REFERENCE_COLOR = "#22313F"                                 # the color of Reference transcripts
MATCH_COLOR = "#52B3D9"                                     # the color of matched isoforms

# Select isoforms of a particular gene
def selectGene(opt, isAnnot, isMatch):
    tranList = list()                                      # list of Transcript objects
    exonList = list()                                      # list of Exon objects
    if isAnnot:
        try:
            getGeneFromAnnotation (opt, tranList, exonList)
        except RuntimeError:
            Console.text = 'Console:\n%s not found in annotation \nfile' % opt.gene
    if isMatch:
        getGeneFromMatches (opt, tranList, exonList)
    return tranList, exonList

def getChromosome(tranList):                                # find out which chromosome does the gene locate on
    chromosome = None
    for tran in tranList:
        if tran.annot is False:                             # only annotation file has this information
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
    geneSource.data = dict(                                             # update the table in visualization
        Gene = allGenes.keys(),
        Isoforms = allGenes.values(),
    )

def updateGene(attrname, old, new):                     # update visualization if the gene changes
    opt.gene = Gene.value.strip()
    matchList = Matches.value.strip().split(',')
    opt.matches = matchList
    p.title = "Transcript of %s" % Gene.value.strip()                           # update the title of plot
    # Reset when updating genes
    blockDict = dict(top=[], bottom=[], left=[], right=[], exon=[],             # all the data when updating plot
                    start=[], end=[], chromosome=[], xs=[], ys=[],
                    hover_fill_color=[], boundary=[])
    sourceDict = dict(name=[], xs=[], ys=[], colors=[], line_alpha=[], width=[], height=[],
                        tran=[], full=[], partial=[], annot=[], QScore=[], start=[], end=[])
    codonDict = dict(x=[], y=[], color=[])
    blockSource.data = blockDict
    source.data = sourceDict
    codonSource.data = codonDict
    # create and display a list of all the genes

    # load the reference transcripts
    Console.text = 'Console:\nReading pickle file...'
    if opt.clusterDict is None:                                                 # if it's the first time to load up pickle file
        try:
            clusterDict = getMatchedIsoforms(getParams(None, matchList, None))
            opt.clusterDict = clusterDict                                       # hold pickle file dictionary in RAM
            howManyIsoforms(clusterDict, matchList)                             # find out how many isoforms for each gene
            isMatch = True                                                      # matched file is present
        except IOError:                                                         # if the file is not found in directory
            Console.text = 'Console:\none of the matched file \n%s is not found' % matchList
            isMatch = False
    else:
        if set(opt.clusterDict.keys()) != set(matchList):            # if the matched files are updated
            try:
                clusterDict = dict()
                clusterDict = getMatchedIsoforms(getParams(None, matchList, None))
                opt.clusterDict = clusterDict
                howManyIsoforms(clusterDict, matchList)
                isMatch = True
            except IOError:
                Console.text = 'Console:\none of the matched file \n%s is not found' % matchList
                isMatch = False
        else:
            isMatch = True

    Console.text = 'Console:\nReading annotation file...'
    if opt.annotations is None:
        try:
            opt.gtf = GTF.value.strip()
            Annotations = getAnnotations(opt)           # get the lists of transcripts in annot file, hold it in RAM for quicker respond when changing the gene
            opt.annotations = Annotations
            isAnnot = True
        except IOError:
            Console.text = 'Console:\nannotations file \n%s is not found' % opt.gtf
            isAnnot = False
    else:
        if opt.gtf != GTF.value.strip():
            try:
                opt.gtf = GTF.value.strip()
                Annotations = getAnnotations(opt)
                opt.annotations = Annotations
                isAnnot = True
            except IOError:
                Console.text = 'Console:\nannotations file \n%s is not found' % opt.gtf
                isAnnot = False
        else:
            isAnnot = True

    global length, colorDF
    tranList, exonList = selectGene(opt, isAnnot, isMatch)                    # select transcripts by gene
    chromosome = getChromosome(tranList)                    # find out which chromosome does the gene locate

    forwardStrand = '-' if opt.flip else '+'
    if exonList[0].strand == forwardStrand:
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = assignBlocks (opt, exonList)              # assign each exon to a block
    else:
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = assignBlocksReverse (opt, exonList)       # assign each exon to a block -- backwards
    for exon in exonList:
        if exon.tran.annot:
            blocks[exon.block].annot = True
    findRegions (tranList)                       # determine regions occupied by each transcript
    tranNames = orderTranscripts (tranList)                 # get the names of transcripts, placed them in the right order
    tranNames = reduceNameLength(tranNames)                      # if the length of name is too long, reduce it

    length = len(tranNames)
    Console.text = 'Console:\nCreating plot...'

    p.plot_height = Height.value*2*(length+4)       # set the height of plot according to the length of transcripts
    p.y_range.factors = tranNames[::-1]             # set the y axis tick to the transcripts names

    Console.text = 'Console:\nGrouping...'
    if Group.value == "on" and isMatch is True:
        colorDF = groupTran(tranList, exonList, 15)          # group the transcripts by similarities
    else:
        colorDF = None
    sourceDict = getExonData(exonList, colorDF)                  # get the data of each isoform that can be directly used to plot

    codonDict = plotStartStop (tranList, blocks)

    codonSource.data = codonDict
    source.data = sourceDict
    # update the data used for plotting boundaries and hover block
    blockDict = getBoundaryData(blocks, chromosome)              # get the data of each block that can be directly used to plot
    blockSource.data = blockDict

    if isAnnot == False:
        Console.text = 'Console:\nSuccess! Annotation\n file is missing.'
    elif isMatch == False:
        Console.text = 'Console:\nSuccess! Match file\n is missing.'
    else:
        Console.text = 'Console:\nSuccess!'

# update plot if full and partial support threshold changes
def updateFP(attrname, old, new):
    sourceDict = source.data
    sourceDict['line_alpha'] = [greaterFP(pos) for pos in range(len(sourceDict['annot']))]
    source.data = sourceDict

# update the number of groups
def updateGroup(attrname, old, new):
    colors = list()
    sourceDict = source.data
    colors = [getColorFromDF(tran, colorDF) for tran in sourceDict['tran']]
    sourceDict['colors'] = colors
    source.data = sourceDict

# update the width of the each isoform
def updateHeightWidth(attrname, old, new):
    p.plot_height = Height.value*2*(length+4)
    p.plot_width = Width.value
    sourceDict = source.data
    sourceDict['height'] = [Height.value for x in range(len(sourceDict['xs']))]
    source.data = sourceDict
    codonDict = codonSource.data
    codonDict['size'] = [Height.value*1.2 for x in range(len(codonDict['x']))]
    codonSource.data = codonDict

# get the data for plotting exons (start, end position for example)
def getExonData(exonList, colorDF):
    sourceDict = dict(name=[], xs=[], ys=[], colors=[], line_alpha=[], width=[], height=[],
                        tran=[], full=[], partial=[], annot=[], QScore=[], start=[], end=[])
    columns = ['name', 'xs', 'ys', 'colors', 'QScore',
                'start', 'end', 'tran', 'full', 'partial', 'annot']
    for myExon in exonList:
        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart
        if colorDF is not None:
            color = getColorFromDF(myExon.tran.name, colorDF)           # find out what group does the exon belongs
        else:                                                           # if the grouping effect is off, paint default color
            if myExon.tran.annot:
                color = REFERENCE_COLOR
            else:
                color = MATCH_COLOR
        xs = (adjStart, adjStart+exonSize)
        ys = (length-(myExon.tran.tranIx), length-(myExon.tran.tranIx))
        values = [myExon.name, xs, ys, color,
                    myExon.QScore, myExon.start, myExon.end,
                    myExon.tran.name, myExon.full, myExon.partial,
                    myExon.tran.annot]
        counter = 0
        while counter < len(columns):
            sourceDict[columns[counter]].append(values[counter])
            counter += 1

    sourceDict['line_alpha'] = [1 for x in range(len(sourceDict['xs']))]
    sourceDict['height'] = [Height.value for x in range(len(sourceDict['xs']))]
    return sourceDict

# get the color for matched isoforms
def getColorFromDF(exonName, colorDF):
    if exonName not in list(colorDF.name):
        color = '#22313F'
    else:
        row = colorDF.loc[colorDF['name'] == exonName]
        colorName = 'color%s' %str(Cluster.value)
        try:
            color = row[colorName]
        except (ValueError, KeyError):
            color = row['color1']
    return color

# find out the position of boundaries
def getBoundaryData(blocks, chromosome):
    blockDict = dict(top=[], bottom=[], left=[], right=[], exon=[],
                    start=[], end=[], chromosome=[], xs=[], ys=[],
                    hover_fill_color=[], boundary=[])
    blockDict['xs'] = [(0, 0)]
    blockDict['ys'] = [(0, length+1)]
    blockDict['left'] = [0]
    columns = ['boundary', 'right', 'left', 'xs', 'start', 'end', 'exon', 'hover_fill_color']
    exonCounter = 1

    for bound in blocks:                                # infomation for the mouse hover effect on blocks
        if bound.annot:
            exon = exonCounter
            color = 'red'
            exonCounter += 1
        else:
            exon = None
            color = 'blue'
        values = [bound.boundary, bound.boundary, bound.boundary,
                        (bound.boundary, bound.boundary), bound.start,
                        bound.end, exon, color]
        counter = 0

        numberOfBlocks = len(blocks)
        blockDict['top'] = [(length+1) for x in range(numberOfBlocks)]
        blockDict['bottom'] = [0 for x in range(numberOfBlocks)]
        blockDict['chromosome'] = [chromosome for x in range(numberOfBlocks)]

        while counter < len(columns):
            blockDict[columns[counter]].append(values[counter])
            counter += 1
    left = blockDict['left']
    blockDict['left'] = left[:-1]
    blockDict['ys'] = [(0, length+1) for x in range(numberOfBlocks+1)]
    return blockDict

# find out which isoform is below full/partial threshold, which isoform is not
def greaterFP(pos):
    alphaVal = Alpha.value
    sourceDict = source.data
    annot = sourceDict['annot']
    if annot[pos] is False:
        if sourceDict['full'][pos] < Full.value or sourceDict['partial'][pos] < Partial.value:
            return alphaVal                 # change alpha values of those below alpha threshold
        else:
            return 1.0
    else:
        return 1.0

# save the transcripts to .fasta file, the function is copied from MatchAnnot
def saveFasta(attrname, old, new):
    Console.text = 'Console:\nSaving...'
    opt.fasta = Save.value.strip()
    tranList = list()
    exonList = list()
    getGeneFromMatches (opt, tranList, exonList)
    opt.fasta = None
    Console.text = 'Console:\nSuccessfully saved'

# create the visualization plot
def createPlot():
    p = Figure(plot_height=300, plot_width=Width.value, title="", y_range=[],
                title_text_font_size=TITLE_FONT_SIZE)
    p.xgrid.grid_line_color = None                              # get rid of the grid in bokeh
    p.ygrid.grid_line_color = None

    quad = p.quad(top="top", bottom="bottom", left="left", right="right", source=blockSource,           # create quad when mouse hover
        fill_color="grey", hover_fill_color="hover_fill_color",
        fill_alpha=0.05, hover_alpha=0.3,
        line_color=None, hover_line_color="white")
    p.multi_line(xs="xs", ys="ys", source=source, color="colors",                            # plot exons
                line_width="height", line_alpha='line_alpha')
    p.multi_line(xs="xs", ys="ys", source=blockSource, color="black",                       # plot boundaries
                line_width=2, line_alpha=0.4, line_dash="dotted")
    p.inverted_triangle(x="x", y="y", color="color", source=codonSource, size='size', alpha=0.5)
    p.add_tools(HoverTool(tooltips=[("chromosome", "@chromosome"),("exon", "@exon"),        # make mouse hover work
                ("start", "@start"), ("end", "@end")], renderers=[quad]))
    return p

def plotStartStop (tranList, blocks):
    '''Add start/stop codons to plot.'''
    codonDict = dict(x=[], y=[], color=[])
    length = len(tranList)
    for tran in tranList:
        if tran.annot:                             # only annotations know about start/stops
            if hasattr(tran, 'startcodon'):
                codonDict['color'].append('green')
                xPos = findCodon(tran.startcodon, blocks)
                codonDict['x'].append(xPos)
                codonDict['y'].append(length - tran.tranIx)
            if hasattr(tran, 'stopcodon'):
                codonDict['color'].append('red')
                xPos = findCodon(tran.stopcodon, blocks)
                codonDict['x'].append(xPos)
                codonDict['y'].append(length - tran.tranIx)
    codonDict['size'] = [Height.value*1.2 for x in range(len(codonDict['x']))]
    return codonDict

def findCodon (posit, blocks):
    '''Add a codon mark to the plot.'''
    for blk in blocks:

        if blk.start <= posit and blk.end >= posit or \
                blk.start >= posit and blk.end <= posit:      # check in both strand directions
            xPos = blk.boundary - abs(blk.end-posit)
    return xPos

# a class containing all the input parameters
class getParams(object):
    def __init__(self, gtf, matches, gene, forMat="standard", omit=None,
                 show=None, howmany=None, nodups=None, minlen=None, maxlen=None, output="exon.png",
                 flip=None, yscale=1.0, details=None, fasta=None, title=None, notes=None, full=None,
                 partial=None, highsupport=None, annotations=None, clusterDict=None):
        self.gtf = gtf
        self.matches = matches
        self.gene = gene
        self.format = forMat                        #format is a keyword
        self.omit = omit
        self.show = show
        self.howmany = howmany
        self.nodups = nodups
        self.minlen = minlen
        self.maxlen = maxlen
        self.output = output
        self.flip = flip
        self.yscale = yscale
        self.details = details
        self.fasta = fasta
        self.title = title
        self.notes = notes
        self.highsupport = highsupport
        self.full = full
        self.partial = partial
        self.annotations = annotations
        self.clusterDict = clusterDict

# create all kinds of widgets
GTF = TextInput(title="Enter the name of annotation file", value="gencode.vM9.annotation.gtf")
Format = TextInput(title="Enter the format of annotation file, standard is gtf", value="standard")
Matches = TextInput(title="Enter the name of pickle files from MatchAnnot,e.g. of multiple files: a.pickle,b.pickle", value="matches.pickle")
Gene = TextInput(title="Select gene to visualize")
Alpha = Slider(title="Alpha value of exons", value=1.0, start=0, end=1.0, step=0.1)
Full = Slider(title="Full support threshold", value=0, start=0, end=30, step=1.0)
Partial = Slider(title="Partial support threshold", value=0, start=0, end=50, step=1.0)
Group = Select(title="Group isoform or not", value="on", options=["on", "off"])
Cluster = Slider(title="The number of groups", value=3, start=1, end=15, step=1.0)
Height = Slider(title="The height of transcripts", value=20, start=5, end=30, step=1)
Width = Slider(title="The width of plot", value=1200, start=400, end=1500, step=50)
Save = TextInput(title="Enter the folder name to data in Fasta", value=None)

# the data used for plotting isoforms, boundaries and gene
matchList = Matches.value.strip().split(',')
opt = getParams(GTF.value.strip(), matchList, None, forMat=Format.value.strip())

blockDict = dict(top=[], bottom=[], left=[], right=[], exon=[],
                start=[], end=[], chromosome=[], xs=[], ys=[],
                hover_fill_color=[], boundary=[])
sourceDict = dict(name=[], xs=[], ys=[], colors=[], line_alpha=[], width=[], height=[],
                    tran=[], full=[], partial=[], annot=[], QScore=[], start=[], end=[])
geneDict = dict(Gene=[], Isoforms=[])
codonDict = dict(x=[], y=[], color=[], size=[])

blockSource = ColumnDataSource(data=blockDict)
source = ColumnDataSource(data=sourceDict)
geneSource = ColumnDataSource(data=geneDict)
codonSource = ColumnDataSource(data=codonDict)
# the console box
Console = PreText(text='Console:\nStart visualize by entering \nannotations, pickle file and gene.\nPress Enter to submit.\n',
                    width=250, height=100)
# the visualization plot
p = createPlot()

# a table of with all the genes in the match files, and how many isoforms in each gene
geneColumns = [
        TableColumn(field="Gene", title="Gene"),
        TableColumn(field="Isoforms", title="Isoforms"),
    ]
geneCountTable = DataTable(source=geneSource, columns=geneColumns, width=200, height=1200)

# a description of all the parameters
paramSource = ColumnDataSource(data=dict(
        Parameter=['annotation', 'format', 'matches', 'gene', 'full', 'alpha', 'partial', 'group', '# of groups', 'fasta'],
        Description=['annotations file, in format specified by --format, Reload page to update' ,
                    'format of annotation file: standard, alt, pickle (default=standard: gtf)',
                    'pickle file from matchAnnot.py, if there are multiple files, seperate them with comma no space, Reload page to update',
                    'gene to plot (required)',
                    'add alpha channel to transcripts without low full/partial support, change full/partial to update',
                    'full support threshold, work together with alpha',
                    'partial support threshold, work together with alpha',
                    'group the trascript on similarity',
                    'assign transcripts into how many groups',
                    'out put the .fasta files to a folder, enter the folder name'],
))

paramColumns = [
        TableColumn(field="Parameter", title="Parameter"),
        TableColumn(field="Description", title="Description"),
    ]
paramTable = DataTable(source=paramSource, columns=paramColumns, width=1200, height=700)

# make changes to the plot when widgets are updated
Gene.on_change('value', updateGene)
Full.on_change('value', updateFP)
Partial.on_change('value', updateFP)
Alpha.on_change('value', updateFP)
Cluster.on_change('value', updateGroup)
Save.on_change('value', saveFasta)
Height.on_change('value', updateHeightWidth)
Width.on_change('value', updateHeightWidth)

# the position of plot and widgets
files = [GTF, Format, Matches]
controls = [Console, Gene, Height, Width, Alpha, Full, Partial, Group, Cluster, Save]
main = VBoxForm(p, HBox(*files, width=1100), paramTable)
inputs = HBox(VBoxForm(*controls), width=250)
curdoc().add_root(HBox(inputs, main, geneCountTable, width=1800))
