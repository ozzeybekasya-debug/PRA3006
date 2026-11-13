# PRA3006 - Programming in Life Sciences

## Description
This is the final product of the practical from the Maastricht Science Programme (Maastricht University). The goal of this project is to interact SPARRQL endpoints with JavaScript and visualize the results graphically using the d3.js library.

## Research Question
To what extent is gene overlap between human metabolic pathways (from Wikipathways) associated with functional similarity among their enzymes encoded by nonshared genes, based on Enzyme Commission (EC) classification? 

## Dataset
Pathway-gene associations were retrieved from Wikidata using SPARQL. Human pathways were selected using their WikiPathways identifiers.

## Data Collection
SPARQL queries are stored in the queries folder. These queries can be used to retrieve data from Wikidata Query Service.

1. pathway_genes.sparql - retrieving pathway associated genes

2. enzyme_ec.sparql - retrieving EC numbers