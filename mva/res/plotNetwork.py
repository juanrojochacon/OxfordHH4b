#! /usr/bin/python
import matplotlib.pyplot as plt
import networkx as nx
import os
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy
import sys

# Init data
sourcefile = 'nn_9X5X3X1_500000-Gen_CE.net'

# Verify paths
if os.path.exists(sourcefile) == False:
  	print "Error: source file" + sourcefile + " not found!"
  	sys.exit()

infile = open(sourcefile, 'rb')
print "Processing " + sourcefile + " ..."

# Read architecture
line = infile.readline()
arch = []
for ls in line.split():
	arch.append(int(ls))
nNodes = sum(arch)

# Read kinematics names
line = infile.readline()
kinematics = []
for ls in line.split():
	kinematics.append(ls)

# Init graph
G=nx.Graph()

weights = []
thresholds = [0]

line = infile.readline()
for line in infile:
	iNode = int(line.split()[0])
	tNode = int(line.split()[1])
	threshold = float(line.split()[2])
	weight = abs(float(line.split()[3]))

	weights.append(weight)
	thresholds.append(threshold)

	# Add the edge
	G.add_edge(iNode, tNode, wgt = weight )
	# Add threshold weight
	G.node[iNode]['bias'] = threshold
	if tNode < arch[0]:
		G.node[tNode]['bias'] = 0

	# Init labels
	G.node[iNode]['label'] = ''
	if tNode < arch[0]:
		G.node[tNode]['label'] = tNode
	else:
		G.node[tNode]['label'] = ''

# Get input weights
sumw = []
for i in xrange(0,arch[0]):
	G.node[i]['sumw'] = 0
	for j in G.neighbors(i):
		G.node[i]['sumw']=G.node[i]['sumw']+G[i][j]['wgt']
	sumw.append(G.node[i]['sumw'])
	print kinematics[i], G.node[i]['sumw']


hot = cm = plt.get_cmap('hot') 
cNorm  = colors.Normalize(vmin=min(weights), vmax=max(weights))
weightMap = cmx.ScalarMappable(norm=cNorm, cmap=hot)


jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=min(thresholds), vmax=max(thresholds))
biasMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


colorList = []
for edge in G.edges(data=True):
  colorVal = weightMap.to_rgba(edge[2]['wgt'])
  colorList.append(colorVal)

biasList = []
for node in G.nodes(data=True):
  	colorVal = biasMap.to_rgba(node[1]['bias'])
  	biasList.append(colorVal)

iLayer = 0
iNode = 1
posLoc = []
labelList = []
for node in G.nodes(data=True):
	inc = 0
	if iLayer > 0:
		inc = (arch[0] - arch[iLayer])/2

	loc = [iLayer, iNode + inc]
	iNode = iNode+1
	if iNode > arch[iLayer]:
		iLayer = iLayer+1
		iNode = 1
	posLoc.append(loc)

	labelList.append(node[1]['label'])



plt.xticks(rotation=-25)
plt.bar(numpy.arange(len(kinematics)), sumw, alpha=0.6)
plt.subplots_adjust(bottom=0.20)

plt.xlabel('Input Variable')
plt.ylabel('Total associated weight')
plt.xticks(numpy.arange(len(kinematics)) + 0.5, kinematics)

plt.savefig("nnweights.pdf")
plt.clf()

nx.draw(G, pos=dict(zip(G.nodes(),posLoc)), labels=dict(zip(G.nodes(),labelList)), edge_color=colorList, node_color=biasList, width=2, with_labels=True)
plt.savefig("nnarch.pdf")
plt.clf()
