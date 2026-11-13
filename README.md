# PRA3006 — Programming in Life Sciences

## Functional Similarity and Gene Overlap Between Human Metabolic Pathways

### Description

This project is part of the Programming in Life Sciences practical at the
Maastricht Science Programme (Maastricht University). The goal is to retrieve
data from the Wikidata SPARQL endpoint, analyze the results in JavaScript,
and visualize them on a webpage using D3.js.

### Research Question

To what extent is gene overlap between **human metabolic pathways** (from Wikipathways) associated with functional similarity among their enzymes encoded by nonshared genes, based on **Enzyme Commission (EC) classification**?

### Data Sources

Pathway–gene associations were retrieved from Wikidata using SPARQL. Human metabolic pathways were selected using their WikiPathways identifiers.
Protein-coding status and EC annotations are queried separately in Wikidata.

### Files

- `queries/pathway_genes.sparql` — Retrieves the genes associated with the
  selected pathways.
- `queries/enzyme_ec.sparql` — Retrieves protein-coding status and EC numbers
  for the selected pathway genes where annotations are available.
- `index.html` — Contains the project webpage, which currently displays
  the title and research question.
- `analyze.js` — The JavaScript file for calculating gene overlap and
  comparing EC annotations among non-shared genes.

### License

This project is licensed under the MIT License. See [LICENSE](LICENSE).