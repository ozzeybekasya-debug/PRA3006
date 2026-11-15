import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

const endpoint = "https://query.wikidata.org/sparql";
const pathwayIds = ["WP534", "WP134", "WP690", "WP78"];
const status = document.querySelector("#trial-status");

// 1. Read a saved query, send it to Wikidata, and return its result rows.
async function runQuery(filename) {
  const queryFile = await fetch(new URL(`./queries/${filename}`, import.meta.url));
  if (!queryFile.ok) {
    throw new Error(`Could not load ${filename} (HTTP ${queryFile.status}).`);
  }

  const url = new URL(endpoint);
  url.searchParams.set("query", await queryFile.text());
  url.searchParams.set("format", "json");

  const response = await fetch(url, {
    headers: { Accept: "application/sparql-results+json" }
  });
  if (!response.ok) {
    throw new Error(`Wikidata request for ${filename} failed (HTTP ${response.status}).`);
  }

  const data = await response.json();
  return data.results.bindings;
}

// 2. Compare sets without counting duplicate genes or EC subclasses.
function intersection(a, b) {
  return new Set([...a].filter(value => b.has(value)));
}

function difference(a, b) {
  return new Set([...a].filter(value => !b.has(value)));
}

function jaccard(a, b) {
  const unionSize = new Set([...a, ...b]).size;
  return unionSize ? intersection(a, b).size / unionSize : null;
}

function ecProfile(genes, ecByGene) {
  const subclasses = new Set();
  let annotatedGenes = 0;

  for (const gene of genes) {
    const codes = ecByGene.get(gene);
    if (!codes?.size) continue;

    annotatedGenes++;
    for (const code of codes) subclasses.add(code);
  }

  return { subclasses, annotatedGenes };
}

// 3. Build the six comparisons.
async function main() {
  const [geneRows, ecRows] = await Promise.all([
    runQuery("pathway_genes.sparql"),
    runQuery("enzyme_ec.sparql")
  ]);

  // Keep a slot for every pilot pathway, even if Wikidata returns no genes.
  const genesByPathway = new Map(pathwayIds.map(id => [id, new Set()]));

  for (const row of geneRows) {
    const genes = genesByPathway.get(row.wpid?.value);
    if (genes && row.gene?.value) genes.add(row.gene.value);
  }

  // Use the gene URI as the key in both queries.
  // For this pilot, keep only complete numeric EC numbers.
  // Example: "5.3.1.1" contributes the subclass "5.3".
  const ecByGene = new Map();

  for (const row of ecRows) {
    const gene = row.gene?.value;
    const ec = row.ecNumber?.value?.trim();
    if (!gene || !ec || !/^\d+(?:\.\d+){3}$/.test(ec)) continue;

    if (!ecByGene.has(gene)) ecByGene.set(gene, new Set());
    ecByGene.get(gene).add(ec.split(".").slice(0, 2).join("."));
  }

  const pairs = [];

  for (let i = 0; i < pathwayIds.length; i++) {
    for (let j = i + 1; j < pathwayIds.length; j++) {
      const [idA, idB] = [pathwayIds[i], pathwayIds[j]];
      const a = genesByPathway.get(idA);
      const b = genesByPathway.get(idB);
      const hasGeneLists = a.size > 0 && b.size > 0;

      const onlyA = difference(a, b);
      const onlyB = difference(b, a);
      const ecA = ecProfile(onlyA, ecByGene);
      const ecB = ecProfile(onlyB, ecByGene);

      pairs.push({
        label: `${idA} vs ${idB}`,
        geneCounts: `${a.size} / ${b.size}`,
        shared: hasGeneLists ? intersection(a, b).size : null,
        overlap: hasGeneLists ? jaccard(a, b) : null,
        coverage: hasGeneLists
          ? `${ecA.annotatedGenes}/${onlyA.size} | ${ecB.annotatedGenes}/${onlyB.size}`
          : "Not assessable",
        similarity: hasGeneLists && ecA.subclasses.size && ecB.subclasses.size
          ? jaccard(ecA.subclasses, ecB.subclasses)
          : null
      });
    }
  }

  drawChart(pairs);
  drawTable(pairs);

  const plotted = pairs.filter(pair => pair.similarity !== null).length;
  status.textContent = `${plotted} of ${pairs.length} pairs have both scores.`
    + " Pairs with missing data are still listed below.";
}

// 4. Plot only pairs with both scores.
function drawChart(pairs) {
  const points = pairs.filter(p => p.overlap !== null && p.similarity !== null);
  if (!points.length) return;

  const x = d3.scaleLinear().domain([0, 1]).range([65, 585]);
  const y = d3.scaleLinear().domain([0, 1]).range([305, 25]);
  const percent = d3.format(".0%");

  const svg = d3.select("#trial-chart").append("svg")
    .attr("viewBox", "0 0 610 355")
    .style("width", "100%")
    .attr("role", "img")
    .attr("aria-label", "Gene overlap versus EC subclass similarity");

  svg.append("g")
    .attr("transform", "translate(0,305)")
    .call(d3.axisBottom(x).ticks(5).tickFormat(percent));

  svg.append("g")
    .attr("transform", "translate(65,0)")
    .call(d3.axisLeft(y).ticks(5).tickFormat(percent));

  svg.append("text")
    .attr("x", 325).attr("y", 345).attr("text-anchor", "middle")
    .text("Gene overlap (Jaccard)");

  svg.append("text")
    .attr("transform", "translate(17,165) rotate(-90)")
    .attr("text-anchor", "middle")
    .text("EC subclass similarity (Jaccard)");

  svg.selectAll("circle")
    .data(points)
    .join("circle")
    .attr("cx", d => x(d.overlap))
    .attr("cy", d => y(d.similarity))
    .attr("r", 6)
    .attr("fill", "#167768")
    .attr("opacity", 0.8)
    .append("title")
    .text(d => `${d.label} | Genes: ${percent(d.overlap)} | EC: ${percent(d.similarity)}`);
}

// 5. Show every pair, even when it can't be plotted.
function drawTable(pairs) {
  const score = value => value === null
    ? "Not assessable"
    : d3.format(".1%")(value);

  const table = d3.select("#trial-table").append("table");
  table.append("caption")
    .text("EC coverage counts annotated non-shared genes out of all non-shared genes on each side.");

  const headings = [
    "Pair", "Genes A / B", "Shared genes", "Gene overlap",
    "EC coverage A | B", "EC similarity"
  ];

  table.append("thead").append("tr").selectAll("th")
    .data(headings).join("th").attr("scope", "col").text(d => d);

  table.append("tbody").selectAll("tr")
    .data(pairs).join("tr").selectAll("td")
    .data(p => [
      p.label, p.geneCounts, p.shared ?? "Not assessable",
      score(p.overlap), p.coverage, score(p.similarity)
    ])
    .join("td").text(d => d);
}

main().catch(error => {
  status.textContent = `Could not load the trial results. ${error.message}`;
  console.error(error);
});