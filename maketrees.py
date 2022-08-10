
import os
os.chdir("/Users/mas4037/Dropbox (Mason Lab)/Hillier")

from ete3 import NCBITaxa
ncbi = NCBITaxa() #download taxa from ncbi
taxonomy = [line.rstrip() for line in open("archaea_phylum_taxa.txt", "r")] #Include query taxonomy
name2taxid = ncbi.get_name_translator(taxonomy)
taxid=name2taxid.values()
taxnames=name2taxid.items()

taxid=sum(taxid, [])
tree = ncbi.get_topology(taxid, intermediate_nodes=True)
ncbi.annotate_tree(tree, taxid_attr="name")

for node in tree.traverse():
    node.name = node.sci_name
tree.add_child(name="Halobacterota") #Added since it is not present in NCBI
tree.add_child(name="Euryarchaeota") #Added since it was included as a node in get_topology
tree.write(outfile="tree.archaea.txt", format=1)


#Extremophile tree
from ete3 import NCBITaxa
ncbi = NCBITaxa() #download taxa from ncbi
taxonomy = [line.rstrip() for line in open("hillier_ex_taxa.txt", "r")] #Include query taxonomy
name2taxid = ncbi.get_name_translator(taxonomy)
taxid=name2taxid.values()
taxnames=name2taxid.items()

taxid=sum(taxid, [])
tree = ncbi.get_topology(taxid, intermediate_nodes=True)
ncbi.annotate_tree(tree, taxid_attr="name")

for node in tree.traverse():
    node.name = node.sci_name
tree.write(outfile="tree_extremophiles.txt", format=1) #Must change spaces with _ before importing in R


