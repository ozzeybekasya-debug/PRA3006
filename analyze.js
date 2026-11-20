let d3;

const endpoint = "https://query.wikidata.org/sparql";
const pathwayIds = [
  "WP534", "WP134", "WP500", "WP78",
  "WP357", "WP143", "WP311", "WP197",
  "WP1495", "WP465", "WP497", "WP100"
];
const ecLevels = [
  { depth: 1, name: "Class" },
  { depth: 2, name: "Subclass" },
  { depth: 3, name: "Sub-subclass" },
  { depth: 4, name: "Full EC number" }
];
const status = document.querySelector("#trial-status");
const heading = document.querySelector("#trial-heading");
if (heading) heading.textContent = "Pathway comparison";

// Load a saved SPARQL query and send it to Wikidata.
async function runQuery(filename) {
  const queryPath = new URL(`./queries/${filename}`, import.meta.url);
  const queryFile = await fetch(queryPath);

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

// Compare sets without double-counting genes or EC entries.
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
  const byLevel = ecLevels.map(() => new Set());
  let annotatedGenes = 0;

  for (const gene of genes) {
    const codes = ecByGene.get(gene);
    if (!codes?.size) continue;

    annotatedGenes++;
    for (const code of codes) {
      const parts = code.split(".");
      // Example: 3.5.4.6 contributes 3, 3.5, 3.5.4, and 3.5.4.6.
      ecLevels.forEach(({ depth }, i) => {
        byLevel[i].add(parts.slice(0, depth).join("."));
      });
    }
  }

  return { byLevel, annotatedGenes };
}

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

  // Keep all selected pathways, even if one returns no gene rows.
  const genesByPathway = new Map(pathwayIds.map(id => [id, new Set()]));

  for (const row of geneRows) {
    const genes = genesByPathway.get(row.wpid?.value);
    if (genes && row.gene?.value) {
      genes.add(row.gene.value);
    }
  }

  // Join annotations to pathway members by gene URI.
  // Use the same complete EC annotations at all four levels.
  const ecByGene = new Map();

  for (const row of ecRows) {
    const gene = row.gene?.value;
    const ec = row.ecNumber?.value?.trim();

    if (!gene || !ec || !/^\d+(?:\.\d+){3}$/.test(ec)) continue;

    if (!ecByGene.has(gene)) {
      ecByGene.set(gene, new Set());
    }

    ecByGene.get(gene).add(ec);
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
        hasGeneLists && ecA.annotatedGenes > 0 && ecB.annotatedGenes > 0;

      pairs.push({
        label: `${idA} vs ${idB}`,
        geneCounts: `${a.size} / ${b.size}`,
        shared: hasGeneLists ? intersection(a, b).size : null,
        overlap: hasGeneLists ? jaccard(a, b) : null,
        coverage: hasGeneLists
          ? `${ecA.annotatedGenes}/${onlyA.size} | ${ecB.annotatedGenes}/${onlyB.size}`
          : "Not assessable",
        similarity: ecLevels.map((_, i) =>
          hasEcProfiles ? jaccard(ecA.byLevel[i], ecB.byLevel[i]) : null
        )
      });
    }
  }

  drawTable(pairs);
  drawCharts(pairs);

  const plotted = pairs.filter(
    pair => pair.overlap !== null && pair.similarity.every(score => score !== null)
  ).length;

  status.textContent = `${plotted} of ${pairs.length} pairs have gene overlap and EC similarity at all four levels. `
    + "All pairs are listed in the table.";

  addDownloadButton();
}

// Use the same pairs and scales in every panel, so the levels are comparable.
function drawCharts(pairs) {
  const points = pairs.filter(
    pair => pair.overlap !== null && pair.similarity.every(score => score !== null)
  );

  const chart = d3.select("#trial-chart");

  if (!points.length) {
    chart.append("p").text(
      "No pairs can be plotted yet. Check the gene counts and EC coverage in the table."
    );
    return;
  }

  const grid = chart.append("div")
    .style("display", "grid")
    .style("grid-template-columns", "repeat(auto-fit, minmax(min(100%, 340px), 1fr))")
    .style("gap", "16px");

  ecLevels.forEach((level, i) => drawChart(grid, points, level, i));
}

function drawChart(grid, points, level, i) {
  const panel = grid.append("section");
  panel.append("h3").text(`EC level ${level.depth} — ${level.name}`);

  const x = d3.scaleLinear().domain([0, 1]).range([65, 585]);
  const y = d3.scaleLinear().domain([0, 1]).range([305, 25]);
  const percent = d3.format(".0%");

  const svg = panel.append("svg")
    .attr("viewBox", "0 0 610 355")
    .style("width", "100%")
    .attr("role", "img")
    .attr("aria-label", `Gene overlap versus EC level ${level.depth} similarity`);

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
    .text("EC similarity (Jaccard)");

  svg.selectAll("circle")
    .data(points)
    .join("circle")
    .attr("cx", pair => x(pair.overlap))
    .attr("cy", pair => y(pair.similarity[i]))
    .attr("r", 6)
    .attr("fill", "#167768")
    .attr("opacity", 0.8)
    .append("title")
    .text(pair =>
      `${pair.label} | Gene overlap: ${percent(pair.overlap)}`
      + ` | EC level ${level.depth} similarity: ${percent(pair.similarity[i])}`
    );
}

// List all pathway pairs, including those that cannot be plotted.
function drawTable(pairs) {
  const percent = d3.format(".1%");
  const score = value => value === null ? "Not assessable" : percent(value);

  const table = d3.select("#trial-table").append("table");

  table.append("caption").text(
    "EC coverage shows how many non-shared genes on each side have a usable EC annotation."
  );

  const headings = [
    "Pair", "Genes A / B", "Shared genes",
    "Gene overlap", "EC coverage A | B",
    ...ecLevels.map(({ depth }) => `EC level ${depth}`)
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
      ...pair.similarity.map(score)
    ])
    .join("td")
    .text(value => value);
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
    link.download = "PRA3006-pilot-results.html";
    document.body.append(link);
    link.click();
    link.remove();

    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });

  document.querySelector("#trial-results").append(button);
}

main().catch(error => {
  status.textContent = `Could not load the trial results. ${error.message}`;
  console.error(error);
});
