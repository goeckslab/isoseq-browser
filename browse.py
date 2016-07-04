from getGene import *
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HoverTool, HBox, VBoxForm, TapTool
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, Select, TextInput, PreText, DataTable, TableColumn, CheckboxGroup
from collections import Counter

PLOT_WIDTH = 1200                                           # the width of plott
TITLE_FONT_SIZE = "25pt"                                    # font size of plot title
COLORS = ["#22313F", '#52B3D9', '#BE90D4', '#446CB3', '#86E2D5', '#F5D76E',
            '#F1A9A0', '#663399', '#87D37C', '#26C281', '#96281B',
            '#4ECDC4', '#F4B350', '#6C7A89', '#C5EFF7', '#EF4836']

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
                    start=[], end=[], chromosome=[], xs=[], ys=[], boundary=[])
    sourceDict = dict(name=[], xs=[], ys=[], colors=[], line_alpha=[], width=[], height=[],
                        tran=[], full=[], partial=[], annot=[], QScore=[], start=[], end=[])
    codonDict = dict(x=[], y=[], color=[])
    blockSource.data = blockDict
    allBlockSource.data = blockDict
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

    global length, colorDF, chromosome, strand
    tranList, exonList = selectGene(opt, isAnnot, isMatch)                    # select transcripts by gene
    chromosome = getChromosome(tranList)                    # find out which chromosome does the gene locate

    strand = exonList[0].strand
    if exonList[0].strand == '+':
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
    if 1 in Group.active and isMatch is True:
        colorDF = groupTran(tranList, exonList, 15)          # group the transcripts by similarities
    else:
        colorDF = None
    sourceDict = getExonData(exonList, colorDF)                  # get the data of each isoform that can be directly used to plot

    codonDict = plotStartStop (tranList, blocks)

    codonSource.data = codonDict
    source.data = sourceDict
    # update the data used for plotting boundaries and hover block
    blockDict, tranDict = getBoundaryData(blocks, chromosome)              # get the data of each block that can be directly used to plot
    blockSource.data = blockDict
    allBlockSource.data = blockDict
    tranSource.data = tranDict
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
    sourceDict = source.data
    if 0 in Group.active:
        sourceDict['colors'] = sourceDict['fileColor']
    else:
        if 1 in Group.active:
            colors = list()
            colors = [getColorFromDF(tran, colorDF) for tran in sourceDict['tran']]
            sourceDict['colors'] = colors
        else:
            colors = list()
            for i in sourceDict['annot']:
                if i is True:
                    colors.append(COLORS[0])
                else:
                    colors.append(COLORS[1])
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
                        tran=[], full=[], partial=[], annot=[], QScore=[], start=[], end=[],
                        fileColor=[])
    columns = ['name', 'xs', 'ys', 'colors', 'QScore',
                'start', 'end', 'tran', 'full', 'partial', 'annot', 'fileColor']
    for myExon in exonList:
        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart

        if 0 in Group.active:
            color = COLORS[myExon.tran.source[0]]
        else:
            if colorDF is not None:
                color = getColorFromDF(myExon.tran.name, colorDF)           # find out what group does the exon belongs
            else:                                                           # if the grouping effect is off, paint default color
                if myExon.tran.annot:
                    color = COLORS[0]
                else:
                    color = COLORS[1]
        xs = (adjStart, adjStart+exonSize)
        ys = (length-(myExon.tran.tranIx), length-(myExon.tran.tranIx))
        values = [myExon.name, xs, ys, color,
                    myExon.QScore, myExon.start, myExon.end,
                    myExon.tran.name, myExon.full, myExon.partial,
                    myExon.tran.annot, COLORS[myExon.tran.source[0]]]

        for ix, col in enumerate(columns):
            sourceDict[columns[ix]].append(values[ix])

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
                    start=[], end=[], chromosome=[], xs=[], ys=[], boundary=[])
    tranDict = dict(top=[], bottom=[], left=[], right=[])
    blockDict['xs'] = [(0, 0)]
    blockDict['ys'] = [(0, length+1)]
    columns = ['boundary', 'left', 'right', 'xs', 'start', 'end', 'exon']
    exonCounter = 1
    numberOfBlocks = len(blocks)
    for bound in blocks:                                # infomation for the mouse hover effect on blocks
        values = [bound.boundary, bound.boundary-bound.start+bound.end, bound.boundary,
                (bound.boundary, bound.boundary), bound.start, bound.end, exonCounter]
        blockDict['top'] = [(length+1) for x in range(numberOfBlocks)]
        blockDict['bottom'] = [0 for x in range(numberOfBlocks)]
        blockDict['chromosome'] = [chromosome for x in range(numberOfBlocks)]
        exonCounter += 1
        counter = 0
        while counter < len(columns):
            blockDict[columns[counter]].append(values[counter])
            counter += 1
    right = blockDict['right']
    # # blockDict['left'] = left[:-1]
    # # blockDict['right'] = left[1:]
    blockDict['ys'] = [(0, length+1) for x in range(numberOfBlocks+1)]

    tranDict['top'] = [x+1.5 for x in range(length)]
    tranDict['bottom'] = [x+0.5 for x in range(length)]
    tranDict['left'] = [0 for x in range(length)]
    tranDict['right'] = [max(right) for x in range(length)]
    return blockDict, tranDict

# find out which isoform is below full/partial threshold, which isoform is not
def greaterFP(pos):
    sourceDict = source.data
    annot = sourceDict['annot']
    if annot[pos] is False:
        if sourceDict['full'][pos] < Full.value or sourceDict['partial'][pos] < Partial.value:
            return 0                 # change alpha values of those below alpha threshold
        else:
            return 1
    else:
        return 1

def selection(attr, old, new):
    if tranSource.selected['1d']['indices'] == []:
        sourceDict = source.data
        sourceDict['line_alpha'] = [1 for y in sourceDict['ys']]
        source.data = sourceDict
        blockSource.data = allBlockSource.data
    else:
        index = tranSource.selected['1d']['indices'][0]
        sourceDict = source.data
        alpha = [selectTran(index, y) for y in sourceDict['ys']]
        sourceDict['line_alpha'] = alpha
        source.data = sourceDict
        blockDict = blockSource.data
        blocks = list()
        exons = list()
        for i, yy in enumerate(sourceDict['ys']):
            if yy[0] == index+1:
                if strand == '+':
                    start = sourceDict['start'][i]
                    end = sourceDict['end'][i]
                else:
                    start = sourceDict['end'][i]
                    end = sourceDict['start'][i]
                boundary = sourceDict['xs'][i][1]
                bound = Block(start, end, boundary)
                blocks.append(bound)
        bd, tr = getBoundaryData(blocks, chromosome)
        blockSource.data = bd

def selectTran(index, y):
    if index+1 == y[0]:
        return 1
    else:
        return 0.2

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
    TOOLS="pan,wheel_zoom,save,reset,tap"
    p = Figure(plot_height=300, plot_width=Width.value, title="", y_range=[],
                title_text_font_size=TITLE_FONT_SIZE, webgl=True, tools=TOOLS)
    p.xgrid.grid_line_color = None                              # get rid of the grid in bokeh
    p.ygrid.grid_line_color = None

    quad = p.quad(top="top", bottom="bottom", left="left", right="right", source=blockSource,           # create quad when mouse hover
        hover_fill_color="red", line_dash="dotted",
        fill_alpha=0, hover_alpha=0.3, line_color='black', hover_line_color="white",
        line_alpha=0.4, nonselection_fill_alpha=0, nonselection_line_alpha=0.4,
        nonselection_line_color='black')
    p.quad(top="top", bottom="bottom", right="right", left="left", source=tranSource,
           fill_alpha=0, line_alpha=0, nonselection_fill_alpha=0, nonselection_line_alpha=0)
    p.multi_line(xs="xs", ys="ys", source=source, color="colors",                            # plot exons
                line_width="height", line_alpha='line_alpha')
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
# gencode.vM9.annotation.gtf
GTF = TextInput(title="Enter the name of annotation file", value="gencode.vM9.annotation.gtf")
Format = TextInput(title="Enter the format of annotation file, standard is gtf", value="standard")
Matches = TextInput(title="Enter the name of pickle files from MatchAnnot,e.g. of multiple files: a.pickle,b.pickle", value="matches.pickle")
Gene = TextInput(title="Select gene to visualize")
Full = Slider(title="Full support threshold", value=0, start=0, end=30, step=1.0)
Partial = Slider(title="Partial support threshold", value=0, start=0, end=50, step=1.0)
Group = CheckboxGroup(labels=["Group by file", "Group by similarity"], active=[1])
Cluster = Slider(title="The number of groups", value=3, start=1, end=15, step=1.0)
Height = Slider(title="The height of transcripts", value=20, start=5, end=30, step=1)
Width = Slider(title="The width of plot", value=1200, start=400, end=1500, step=50)
Save = TextInput(title="Enter the folder name to data in Fasta", value=None)

# the data used for plotting isoforms, boundaries and gene
matchList = Matches.value.strip().split(',')
opt = getParams(GTF.value.strip(), matchList, None, forMat=Format.value.strip())

blockDict = dict(top=[], bottom=[], left=[], right=[], exon=[],
                start=[], end=[], chromosome=[], xs=[], ys=[], boundary=[])
tranDict = dict(top=[], bottom=[], left=[], right=[])
sourceDict = dict(exon=[], name=[], xs=[], ys=[], colors=[], line_alpha=[], width=[], height=[],
                    tran=[], full=[], partial=[], annot=[], QScore=[], start=[], end=[],
                    fileColor=[])
geneDict = dict(Gene=[], Isoforms=[])
codonDict = dict(x=[], y=[], color=[], size=[])

blockSource = ColumnDataSource(data=blockDict)
allBlockSource = ColumnDataSource(data=blockDict)
tranSource = ColumnDataSource(data=tranDict)
source = ColumnDataSource(data=sourceDict)
geneSource = ColumnDataSource(data=geneDict)
codonSource = ColumnDataSource(data=codonDict)
# the console box
Console = PreText(text='Console:\nStart visualize by entering \nannotations, pickle file and gene.\nPress Enter to submit.\n',
                    width=250, height=70)
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
Cluster.on_change('value', updateGroup)
Group.on_change('active', updateGroup)
Save.on_change('value', saveFasta)
Height.on_change('value', updateHeightWidth)
Width.on_change('value', updateHeightWidth)
tranSource.on_change('selected', selection)

# the position of plot and widgets
files = [GTF, Format, Matches]
controls = [Console, Gene, Height, Width, Full, Partial, Group, Cluster, Save]
main = VBoxForm(p, HBox(*files, width=1100), paramTable)
inputs = HBox(VBoxForm(*controls), width=250)
curdoc().add_root(HBox(inputs, main, geneCountTable, width=1800))
