# Reactome2Network
Parsing out the reactome information to get it in a good format for cpAOP

This is code used to go from the Reactome data to get the information needed to create the cpAOP network edges and node labels

Data: Reactome pathway information (V53) obtained from http://www.reactome.org/pages/download-data/ (accessed 7/2015)

Code has Two parts

ReactomeClassv2.R -- this takes the "UniProt2Reactome", "ReactomePathwaysRelation" and the  "ReactomePathways" files from reactome and generates parent-child relationships in a table format and brings the uniprot ids in with the pathway annotations

ReactomeXSpecies.R -- Takes the "ReactomePathways" file and reshapes the table to allow for comparing the reactome IDs for the different species for the same pathway
