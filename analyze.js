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

if (heading) {
  heading.textContent = "Pathway comparison";
}


// ======================================================
// 1. RETRIEVE DATA FROM WIKIDATA
// ======================================================

// Load a saved SPARQL query and send it to Wikidata.
async function runQuery(filename) {
  const queryPath = new URL(
    `./queries/${filename}`,
    import.meta.url
  );

  const queryFile = await fetch(queryPath);

  if (!queryFile.ok) {
    throw new Error(
      `Could not load queries/${filename} (HTTP ${queryFile.status}).`
    );
  }

  const url = new URL(endpoint);

  url.searchParams.set("query", await queryFile.text());
  url.searchParams.set("format", "json");

  try {
    const response = await fetch(url, {
      headers: {
        Accept: "application/sparql-results+json"
      },
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

    throw new Error(
      `Wikidata request for ${filename} failed. ${reason}`
    );
  }
}


// ======================================================
// 2. SET OPERATIONS AND JACCARD SIMILARITY
// ======================================================

// Find entries shared by two sets.
function intersection(a, b) {
  return new Set(
    [...a].filter(value => b.has(value))
  );
}

// Find entries in a that are absent from b.
function difference(a, b) {
  return new Set(
    [...a].filter(value => !b.has(value))
  );
}

// Jaccard similarity = intersection / union.
function jaccard(a, b) {
  const unionSize = new Set([...a, ...b]).size;

  return unionSize > 0
    ? intersection(a, b).size / unionSize
    : null;
}


// ======================================================
// 3. BUILD EC PROFILES
// ======================================================

function ecProfile(genes, ecByGene) {
  const byLevel = ecLevels.map(() => new Set());

  let annotatedGenes = 0;

  for (const gene of genes) {
    const codes = ecByGene.get(gene);

    if (!codes?.size) continue;

    annotatedGenes++;

    for (const code of codes) {
      const parts = code.split(".");

      // Example: 3.5.4.6 contributes:
      // Level 1: 3
      // Level 2: 3.5
      // Level 3: 3.5.4
      // Level 4: 3.5.4.6

      ecLevels.forEach(({ depth }, i) => {
        byLevel[i].add(
          parts.slice(0, depth).join(".")
        );
      });
    }
  }

  return {
    byLevel,
    annotatedGenes
  };
}


// ======================================================
// 4. MAIN ANALYSIS
// ======================================================

async function main() {
  status.textContent = "Loading D3…";

  try {
    d3 = await import(
      "https://cdn.jsdelivr.net/npm/d3@7/+esm"
    );
  } catch (error) {
    throw new Error(
      `Could not load D3 from jsDelivr. ${error.message}`
    );
  }

  status.textContent =
    "Loading pathway genes and EC numbers from Wikidata…";

  const [geneRows, ecRows] = await Promise.all([
    runQuery("pathway_genes.sparql"),
    runQuery("enzyme_ec.sparql")
  ]);

  // Keep all 12 selected pathways.
  const genesByPathway = new Map(
    pathwayIds.map(id => [id, new Set()])
  );

  for (const row of geneRows) {
    const genes = genesByPathway.get(
      row.wpid?.value
    );

    if (genes && row.gene?.value) {
      genes.add(row.gene.value);
    }
  }

  // Join EC annotations to pathway genes by gene URI.
  // Only complete four-part EC numbers are accepted.
  const ecByGene = new Map();

  for (const row of ecRows) {
    const gene = row.gene?.value;
    const ec = row.ecNumber?.value?.trim();

    if (
      !gene ||
      !ec ||
      !/^\d+(?:\.\d+){3}$/.test(ec)
    ) {
      continue;
    }

    if (!ecByGene.has(gene)) {
      ecByGene.set(gene, new Set());
    }

    ecByGene.get(gene).add(ec);
  }

  const pairs = [];

  // Compare every unique pair of pathways.
  for (let i = 0; i < pathwayIds.length; i++) {
    for (let j = i + 1; j < pathwayIds.length; j++) {

      const idA = pathwayIds[i];
      const idB = pathwayIds[j];

      const a = genesByPathway.get(idA);
      const b = genesByPathway.get(idB);

      const hasGeneLists =
        a.size > 0 && b.size > 0;

      // Remove shared genes before comparing EC profiles.
      const onlyA = difference(a, b);
      const onlyB = difference(b, a);

      const ecA = ecProfile(onlyA, ecByGene);
      const ecB = ecProfile(onlyB, ecByGene);

      const hasEcProfiles =
        hasGeneLists &&
        ecA.annotatedGenes > 0 &&
        ecB.annotatedGenes > 0;

      pairs.push({
        label: `${idA} vs ${idB}`,

        geneCounts: `${a.size} / ${b.size}`,

        shared: hasGeneLists
          ? intersection(a, b).size
          : null,

        overlap: hasGeneLists
          ? jaccard(a, b)
          : null,

        coverage: hasGeneLists
          ? `${ecA.annotatedGenes}/${onlyA.size} | ${ecB.annotatedGenes}/${onlyB.size}`
          : "Not assessable",

        similarity: ecLevels.map((_, i) =>
          hasEcProfiles
            ? jaccard(ecA.byLevel[i], ecB.byLevel[i])
            : null
        )
      });
    }
  }

  // Generate the results table and scatterplots.
  drawTable(pairs);
  drawCharts(pairs);

  const plotted = pairs.filter(
    pair =>
      pair.overlap !== null &&
      pair.similarity.every(score => score !== null)
  ).length;

  status.textContent =
    `${plotted} of ${pairs.length} pairs have gene overlap ` +
    "and EC similarity at all four levels. " +
    "All pairs are listed in the table.";

  addDownloadButton();
}


// ======================================================
// 5. PEARSON CORRELATION AND LINEAR REGRESSION
// ======================================================

// Pearson's r measures linear association.
//
// The ordinary least-squares regression line is:
//
// y = intercept + slope * x
//
// For simple linear regression with an intercept:
//
// R² = r²

function linearStatistics(points, levelIndex) {
  const n = points.length;

  if (n < 3) {
    return {
      n,
      r: null,
      rSquared: null,
      slope: null,
      intercept: null
    };
  }

  const xMean = d3.mean(
    points,
    pair => pair.overlap
  );

  const yMean = d3.mean(
    points,
    pair => pair.similarity[levelIndex]
  );

  let sxx = 0;
  let syy = 0;
  let sxy = 0;

  for (const pair of points) {
    const dx = pair.overlap - xMean;

    const dy =
      pair.similarity[levelIndex] - yMean;

    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }

  // Pearson correlation is undefined if all x or
  // all y values are identical.
  const epsilon = 1e-14;

  if (sxx <= epsilon || syy <= epsilon) {
    return {
      n,
      r: null,
      rSquared: null,
      slope: null,
      intercept: null
    };
  }

  // Calculate linear regression coefficients.
  const slope = sxy / sxx;

  const intercept =
    yMean - slope * xMean;

  // Calculate Pearson's correlation coefficient.
  const r = Math.max(
    -1,
    Math.min(
      1,
      sxy / Math.sqrt(sxx * syy)
    )
  );

  // Coefficient of determination.
  const rSquared = r * r;

  return {
    n,
    r,
    rSquared,
    slope,
    intercept
  };
}


// ======================================================
// 6. CREATE THE FOUR SCATTERPLOTS
// ======================================================

// Use the same assessable pairs and axis scales
// for all four EC classification levels.

function drawCharts(pairs) {
  const points = pairs.filter(
    pair =>
      pair.overlap !== null &&
      pair.similarity.every(score => score !== null)
  );

  const chart = d3.select("#trial-chart");

  chart.selectAll("*").remove();

  if (!points.length) {
    chart.append("p").text(
      "No pairs can be plotted yet. " +
      "Check the gene counts and EC coverage in the table."
    );

    return;
  }

  const grid = chart.append("div")
    .style("display", "grid")
    .style(
      "grid-template-columns",
      "repeat(auto-fit, minmax(min(100%, 340px), 1fr))"
    )
    .style("gap", "16px");

  ecLevels.forEach((level, i) => {
    drawChart(grid, points, level, i);
  });
}


// ======================================================
// 7. DRAW EACH SCATTERPLOT WITH REGRESSION LINE
// ======================================================

function drawChart(grid, points, level, i) {

  const panel = grid.append("section")
    .style("background-color", "white")
    .style("border", "1px solid #dce8e2")
    .style("border-radius", "12px")
    .style("padding", "14px");

  // Panel heading.
  panel.append("h3")
    .style("margin", "0 0 8px")
    .style("color", "#163b3d")
    .text(
      `EC level ${level.depth} — ${level.name}`
    );

  // Calculate Pearson r, R² and regression coefficients.
  const stats = linearStatistics(points, i);

  const fmt = d3.format(".3f");

  // Display correlation statistics above the chart.
  panel.append("p")
    .style("font-size", "0.9rem")
    .style("margin", "0 0 8px")
    .text(
      stats.r === null
        ? `Pearson r: not defined · R²: not defined · n = ${stats.n}`
        : `Pearson r = ${fmt(stats.r)} · R² = ${fmt(stats.rSquared)} · n = ${stats.n}`
    );

  // Both axes use the same 0–1 Jaccard scale.
  const x = d3.scaleLinear()
    .domain([0, 1])
    .range([65, 585]);

  const y = d3.scaleLinear()
    .domain([0, 1])
    .range([305, 25]);

  const percent = d3.format(".0%");

  const svg = panel.append("svg")
    .attr("viewBox", "0 0 610 355")
    .style("width", "100%")
    .style("display", "block")
    .attr("role", "img")
    .attr(
      "aria-label",
      `Gene overlap and EC level ${level.depth} similarity with linear regression`
    );

  // X-axis.
  svg.append("g")
    .attr("transform", "translate(0,305)")
    .call(
      d3.axisBottom(x)
        .ticks(5)
        .tickFormat(percent)
    );

  // Y-axis.
  svg.append("g")
    .attr("transform", "translate(65,0)")
    .call(
      d3.axisLeft(y)
        .ticks(5)
        .tickFormat(percent)
    );

  // X-axis title.
  svg.append("text")
    .attr("x", 325)
    .attr("y", 345)
    .attr("text-anchor", "middle")
    .attr("fill", "#273841")
    .text("Gene overlap (Jaccard index)");

  // Y-axis title.
  svg.append("text")
    .attr(
      "transform",
      "translate(17,165) rotate(-90)"
    )
    .attr("text-anchor", "middle")
    .attr("fill", "#273841")
    .text("EC similarity (Jaccard index)");

  // Clip regression lines to the plotting area.
  const clipId = `trial-plot-clip-${i}`;

  svg.append("defs")
    .append("clipPath")
    .attr("id", clipId)
    .append("rect")
    .attr("x", 65)
    .attr("y", 25)
    .attr("width", 520)
    .attr("height", 280);

  const plot = svg.append("g")
    .attr(
      "clip-path",
      `url(#${clipId})`
    );

  // --------------------------------------------------
  // ORANGE DASHED LINE: LINEAR REGRESSION
  // --------------------------------------------------

  if (stats.slope !== null) {

    // Draw only within the observed x range.
    const xMin = d3.min(
      points,
      pair => pair.overlap
    );

    const xMax = d3.max(
      points,
      pair => pair.overlap
    );

    const fittedY = value =>
      stats.intercept + stats.slope * value;

    plot.append("line")
      .attr("x1", x(xMin))
      .attr("y1", y(fittedY(xMin)))
      .attr("x2", x(xMax))
      .attr("y2", y(fittedY(xMax)))
      .attr("stroke", "#c46b43")
      .attr("stroke-width", 3)
      .attr("stroke-dasharray", "9 6");
  }

  // --------------------------------------------------
  // GREEN DOTS: PATHWAY PAIRS
  // --------------------------------------------------

  plot.selectAll("circle")
    .data(points)
    .join("circle")
    .attr(
      "cx",
      pair => x(pair.overlap)
    )
    .attr(
      "cy",
      pair => y(pair.similarity[i])
    )
    .attr("r", 6)
    .attr("fill", "#167768")
    .attr("opacity", 0.8)
    .append("title")
    .text(
      pair =>
        `${pair.label} | ` +
        `Gene overlap: ${percent(pair.overlap)} | ` +
        `EC level ${level.depth} similarity: ` +
        percent(pair.similarity[i])
    );

  // ==================================================
  // 8. CHART LEGEND
  // ==================================================

  const legend = panel.append("div")
    .style("display", "flex")
    .style("flex-wrap", "wrap")
    .style("align-items", "center")
    .style("gap", "8px 18px")
    .style("font-size", "0.8rem")
    .style("color", "#526269")
    .style("margin-top", "8px");

  // Green dot legend.
  const dots = legend.append("span")
    .style("display", "inline-flex")
    .style("align-items", "center")
    .style("gap", "7px");

  dots.append("span")
    .style("display", "inline-block")
    .style("width", "10px")
    .style("height", "10px")
    .style("border-radius", "50%")
    .style("background-color", "#167768");

  dots.append("span")
    .text("Pathway pair (Jaccard indices)");

  // Orange dashed regression line legend.
  if (stats.slope !== null) {

    const line = legend.append("span")
      .style("display", "inline-flex")
      .style("align-items", "center")
      .style("gap", "7px");

    line.append("span")
      .style("display", "inline-block")
      .style("width", "28px")
      .style(
        "border-top",
        "3px dashed #c46b43"
      );

    line.append("span")
      .text("Linear regression (best fit)");
  }

  // Explain the Jaccard similarity index.
  panel.append("p")
    .style("font-size", "0.8rem")
    .style("color", "#526269")
    .style("margin", "8px 0 0")
    .text(
      "Jaccard index: 0 = no shared entries; " +
      "1 = identical sets."
    );
}


// ======================================================
// 9. RESULTS TABLE
// ======================================================

// List all pathway pairs, including those that
// cannot be plotted.

function drawTable(pairs) {

  const percent = d3.format(".1%");

  const score = value =>
    value === null
      ? "Not assessable"
      : percent(value);

  const container = d3.select("#trial-table");

  container.selectAll("*").remove();

  const table = container.append("table");

  table.append("caption").text(
    "EC coverage shows how many non-shared genes " +
    "on each side have a usable EC annotation."
  );

  const headings = [
    "Pair",
    "Genes A / B",
    "Shared genes",
    "Gene overlap",
    "EC coverage A | B",
    ...ecLevels.map(
      ({ depth }) => `EC level ${depth}`
    )
  ];

  // Table header.
  table.append("thead")
    .append("tr")
    .selectAll("th")
    .data(headings)
    .join("th")
    .attr("scope", "col")
    .text(heading => heading);

  // Table body.
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


// ======================================================
// 10. DOWNLOAD RESULTS AS HTML
// ======================================================

// Save the rendered webpage, including SVG plots
// and the results table.

function addDownloadButton() {

  const button = document.createElement("button");

  button.id = "download-html";
  button.type = "button";
  button.textContent = "Download HTML";

  button.style.marginTop = "16px";

  button.addEventListener("click", () => {

    const page =
      document.documentElement.cloneNode(true);

    // The copied page already contains the results.
    // Therefore, scripts are no longer needed.
    page.querySelectorAll(
      "script, #download-html"
    ).forEach(element => element.remove());

    const html =
      "<!DOCTYPE html>\n" + page.outerHTML;

    const file = new Blob(
      [html],
      {
        type: "text/html;charset=utf-8"
      }
    );

    const url = URL.createObjectURL(file);

    const link = document.createElement("a");

    link.href = url;
    link.download = "PRA3006-pilot-results.html";

    document.body.append(link);

    link.click();
    link.remove();

    setTimeout(
      () => URL.revokeObjectURL(url),
      1000
    );
  });

  document.querySelector(
    "#trial-results"
  ).append(button);
}


// ======================================================
// 11. RUN ANALYSIS
// ======================================================

main().catch(error => {

  status.textContent =
    `Could not load the trial results. ${error.message}`;

  console.error(error);
});
