# PRA3006 - Programming in Life Sciences

## Research Question
To what extent is gene overlap between human metabolic pathways (from Wikipathways) associated with functional similarity among their enzymes encoded by nonshared genes, based on EC classification? 

## Data Sources
- Wikidata
- WikiPathways
- ExPASy ENZYME

## Data Collection
SPARQL queries are stored in the queries folder. These queries can be used to retrieve data from Wikidata Query Service.

1. pathways.sparql - retrieving human pathways

2. pathway_genes.sparql - retrieving pathway associated genes

3. enzyme_ec.sparql - retrieving EC numbers