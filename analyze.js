let d3;

const endpoint = "https://query.wikidata.org/sparql";
const pathwayIds = ["WP534", "WP134", "WP690", "WP78"];
const status = document.querySelector("#trial-status");

<<<<<<< HEAD
// Load a saved SPARQL query and send it to Wikidata.
async function runQuery(filename) {
  const queryPath = new URL(`./queries/${filename}`, import.meta.url);
  const queryFile = await fetch(queryPath);

=======
// 1. Load a saved SPARQL query and send it to Wikidata.
async function runQuery(filename) {
  const queryPath = new URL(`./queries/${filename}`, import.meta.url);

  const queryFile = await fetch(queryPath);
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
  if (!queryFile.ok) {
    throw new Error(`Could not load queries/${filename} (HTTP ${queryFile.status}).`);
  }

  const url = new URL(endpoint);
  url.searchParams.set("query", await queryFile.text());
  url.searchParams.set("format", "json");

  try {
    const response = await fetch(url, {
      headers: { Accept: "application/sparql-results+json" },
      signal: AbortSignal.timeout(60000)
    });

    if (!response.ok) {
      throw new Error(`HTTP ${response.status}`);
    }

    const data = await response.json();
    return data.results.bindings;
  } catch (error) {
    const reason = error.name === "TimeoutError"
      ? "No response within 60 seconds. Try reloading the page."
      : error.message;

    throw new Error(`Wikidata request for ${filename} failed. ${reason}`);
  }
}

// Compare sets without double-counting genes or EC subclasses.
function intersection(a, b) {
  return new Set([...a].filter(value => b.has(value)));
}

function difference(a, b) {
  return new Set([...a].filter(value => !b.has(value)));
}

function jaccard(a, b) {
  const unionSize = new Set([...a, ...b]).size;
  return unionSize > 0 ? intersection(a, b).size / unionSize : null;
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

<<<<<<< HEAD
=======
// 3. Load D3 and the data, then build the six pairwise comparisons.
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
async function main() {
  status.textContent = "Loading D3…";

  try {
    d3 = await import("https://cdn.jsdelivr.net/npm/d3@7/+esm");
  } catch (error) {
    throw new Error(`Could not load D3 from jsDelivr. ${error.message}`);
  }

  status.textContent = "Loading pathway genes and EC numbers from Wikidata…";

  const [geneRows, ecRows] = await Promise.all([
    runQuery("pathway_genes.sparql"),
    runQuery("enzyme_ec.sparql")
  ]);

<<<<<<< HEAD
  // Keep all four pathways, even if a pathway returns no gene rows.
  const genesByPathway = new Map(pathwayIds.map(id => [id, new Set()]));
=======
  // Keep an entry for each pathway, even if no genes were retrieved.
  const genesByPathway = new Map(
    pathwayIds.map(id => [id, new Set()])
  );
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558

  for (const row of geneRows) {
    const genes = genesByPathway.get(row.wpid?.value);
    if (genes && row.gene?.value) {
      genes.add(row.gene.value);
    }
  }

<<<<<<< HEAD
  // Join annotations to pathway members by gene URI.
  // For this pilot, keep complete numeric EC numbers only.
  // EC 5.3.1.1 contributes the subclass 5.3.
=======
  // Match the two query results by gene URI.
  // For this pilot, use complete numeric EC numbers only.
  // Example: EC 5.3.1.1 contributes the subclass 5.3.
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
  const ecByGene = new Map();

  for (const row of ecRows) {
    const gene = row.gene?.value;
    const ec = row.ecNumber?.value?.trim();

    if (!gene || !ec || !/^\d+(?:\.\d+){3}$/.test(ec)) continue;

    if (!ecByGene.has(gene)) {
      ecByGene.set(gene, new Set());
    }

    ecByGene.get(gene).add(ec.split(".").slice(0, 2).join("."));
  }

  const pairs = [];

  for (let i = 0; i < pathwayIds.length; i++) {
    for (let j = i + 1; j < pathwayIds.length; j++) {
      const idA = pathwayIds[i];
      const idB = pathwayIds[j];
      const a = genesByPathway.get(idA);
      const b = genesByPathway.get(idB);

      const hasGeneLists = a.size > 0 && b.size > 0;
      const onlyA = difference(a, b);
      const onlyB = difference(b, a);
      const ecA = ecProfile(onlyA, ecByGene);
      const ecB = ecProfile(onlyB, ecByGene);
      const hasEcProfiles =
        hasGeneLists && ecA.subclasses.size > 0 && ecB.subclasses.size > 0;

      const hasEcProfiles =
        hasGeneLists && ecA.subclasses.size > 0 && ecB.subclasses.size > 0;

      pairs.push({
        label: `${idA} vs ${idB}`,
        geneCounts: `${a.size} / ${b.size}`,
        shared: hasGeneLists ? intersection(a, b).size : null,
        overlap: hasGeneLists ? jaccard(a, b) : null,
        coverage: hasGeneLists
          ? `${ecA.annotatedGenes}/${onlyA.size} | ${ecB.annotatedGenes}/${onlyB.size}`
          : "Not assessable",
        similarity: hasEcProfiles
          ? jaccard(ecA.subclasses, ecB.subclasses)
          : null
      });
    }
  }

  drawTable(pairs);
  drawChart(pairs);

  const plotted = pairs.filter(
    pair => pair.overlap !== null && pair.similarity !== null
  ).length;

  status.textContent = `${plotted} of ${pairs.length} pairs have both scores. `
    + "All pairs are listed in the table.";
<<<<<<< HEAD

  addDownloadButton();
}

// Plot only pathway pairs that have both similarity scores.
=======
}

// 4. Plot only pairs that have both a gene score and an EC score.
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
function drawChart(pairs) {
  const points = pairs.filter(
    pair => pair.overlap !== null && pair.similarity !== null
  );

  const chart = d3.select("#trial-chart");

<<<<<<< HEAD
  if (!points.length) {
    chart.append("p").text(
      "No pairs can be plotted yet. Check the gene counts and EC coverage in the table."
=======
  if (points.length === 0) {
    chart.append("p").text(
      "No pairs can be plotted yet. Check the EC coverage in the table below."
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
    );
    return;
  }

  const x = d3.scaleLinear().domain([0, 1]).range([65, 585]);
  const y = d3.scaleLinear().domain([0, 1]).range([305, 25]);
  const percent = d3.format(".0%");

  const svg = chart.append("svg")
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
    .attr("x", 325)
    .attr("y", 345)
    .attr("text-anchor", "middle")
    .text("Gene overlap (Jaccard)");

  svg.append("text")
    .attr("transform", "translate(17,165) rotate(-90)")
    .attr("text-anchor", "middle")
    .text("EC subclass similarity (Jaccard)");

  svg.selectAll("circle")
    .data(points)
    .join("circle")
    .attr("cx", pair => x(pair.overlap))
    .attr("cy", pair => y(pair.similarity))
    .attr("r", 6)
    .attr("fill", "#167768")
    .attr("opacity", 0.8)
    .append("title")
    .text(pair =>
      `${pair.label} | Gene overlap: ${percent(pair.overlap)}`
      + ` | EC similarity: ${percent(pair.similarity)}`
    );
}

<<<<<<< HEAD
// List all six pairs, including those that cannot be plotted.
=======
// 5. Show all six pairs, including those that could not be plotted.
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
function drawTable(pairs) {
  const percent = d3.format(".1%");
  const score = value => value === null ? "Not assessable" : percent(value);

  const table = d3.select("#trial-table").append("table");

  table.append("caption").text(
    "EC coverage shows how many non-shared genes on each side have a usable EC annotation."
  );

  const headings = [
<<<<<<< HEAD
    "Pair", "Genes A / B", "Shared genes",
    "Gene overlap", "EC coverage A | B", "EC similarity"
=======
    "Pair",
    "Genes A / B",
    "Shared genes",
    "Gene overlap",
    "EC coverage A | B",
    "EC similarity"
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
  ];

  table.append("thead").append("tr")
    .selectAll("th")
    .data(headings)
    .join("th")
    .attr("scope", "col")
    .text(heading => heading);

  table.append("tbody")
    .selectAll("tr")
    .data(pairs)
    .join("tr")
    .selectAll("td")
    .data(pair => [
      pair.label,
      pair.geneCounts,
      pair.shared ?? "Not assessable",
      score(pair.overlap),
      pair.coverage,
      score(pair.similarity)
    ])
    .join("td")
    .text(value => value);
<<<<<<< HEAD
}

// Save the rendered page, including the SVG chart and results table.
function addDownloadButton() {
  const button = document.createElement("button");
  button.id = "download-html";
  button.type = "button";
  button.textContent = "Download HTML";
  button.style.marginTop = "16px";

  button.addEventListener("click", () => {
    const page = document.documentElement.cloneNode(true);

    // The copy already contains the results, so it doesn't need scripts.
    page.querySelectorAll("script, #download-html")
      .forEach(element => element.remove());

    const html = "<!DOCTYPE html>\n" + page.outerHTML;
    const file = new Blob([html], { type: "text/html;charset=utf-8" });
    const url = URL.createObjectURL(file);

    const link = document.createElement("a");
    link.href = url;
    link.download = "index.html";
    document.body.append(link);
    link.click();
    link.remove();

    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });

  document.querySelector("#trial-results").append(button);
=======
>>>>>>> 316fa13ed5eac3f0ff69eaeaa19b35c78d5e6558
}

main().catch(error => {
  status.textContent = `Could not load the trial results. ${error.message}`;
  console.error(error);
});
