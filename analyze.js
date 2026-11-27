// Initial Setup
// Declaring a variable for the D3.js library to be loaded to in the future.
let d3;

// Defining a constant for the Wikidata endpoint.
const endpoint = "https://query.wikidata.org/sparql"

// Creating an array of four objects for the EC classification levels.
const ecLevels = [
  {depth: 1, name: "Main Enzyme Class"},
  {depth: 2, name: "Subclass"},
  {depth: 3, name: "Sub-Subclass"},
  {depth: 4, name: "Enzyme Identifier"}
]

// Selecting the HTML status element which communicates what the program is currently doing.
const statusTrial = document.querySelector("#trial-status");

// Updating the analysis status if the HTML element exists.
function updateStatus(message) {
  if (statusTrial) {
    statusTrial.textContent = message;
  }
}

// Retrieving Data from Wikidata
// Turning the enpoint string to an URL.
const endpointUrl = new URL(endpoint)

async function runQuery(filename) {
  // Constructing the absolute URL of the SPARQL query file.
  const queryPath = new URL(
    `./queries/${filename}`,
    import.meta.url
  );

  // Retrieving the SPARQL query file.
  const queryFile = await fetch(queryPath);

  // Throwing an error if the query file cannot be loaded.
  if (!queryFile.ok) {
    throw new Error(
      `Could not load queries/${filename} (HTTP ${queryFile.status}).`
    );
  }

  // Attaching the SPARQL query and request JSON formatted results.
  endpointUrl.searchParams.set("query", await queryFile.text());
  endpointUrl.searchParams.set("format", "json");

  try {
    // Sending the query to Wikidata with a JSON response preference and a 100-second timeout.
    const response = await fetch(endpointUrl, {
      headers: {
        Accept: "application/sparql-results+json"
      },
      signal: AbortSignal.timeout(100000)
    });

    // Throwing an error if Wikidata returns an unsuccessful HTTP response.
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}`);
    }

    // Converting the JSON response into a JavaScript object.
    const data = await response.json();

    // Returning the array containing the SPARQL result rows.
    return data.results.bindings;

  }
  
  catch (error) {
    // Distinguishing timeout errors from other request errors.
    const reason = error.name === "TimeoutError"
      ? "No response within 100 seconds. Try reloading the page."
      : error.message;

    // Reporting the failed query and its error message.
    throw new Error(
      `Wikidata request for ${filename} failed. ${reason}`
    );
  }
}

// Jaccard Index
// Returning a Set containing the elements shared by both Sets.
function intersection(a,b) {
  const setA = new Set(a);
  const setB = new Set(b);
  const shared = [...setA].filter(value => setB.has(value));
  return new Set(shared);
}

// Returning a Set containing elements in a that are absent from b.
function difference(a, b) {
  const uniquetoA = new Set(a).difference(new Set(b));
  return uniquetoA;
}

// Calculating the Jaccard similarity between two Sets.
function jaccard(a, b) {
  const union = new Set(a).union(new Set(b)).size;

  if (union === 0) {
    return null;
  }

  const similarity = intersection(a, b).size / union;
  return similarity;
}

// EC Profiles at 4 Levels
// Constructing EC profiles at all four classification levels for a given Set of genes.
function ecProfile(genes, ecByGene) {

  // Creating four independent empty Sets to store unique EC classifications.
  const byLevel = ecLevels.map(() => new Set());

  // Initializing a counter for genes with usable EC annotations.
  let annotatedGenes = 0;

  // Iterating over every gene in the provided Set.
  for (const gene of genes) {

    // Retrieving the EC numbers associated with the current gene.
    const codes = ecByGene.get(gene);

    // Skipping the gene if it has no usable EC annotations.
    if (!codes?.size) continue;

    // Counting the gene as having usable EC annotations.
    annotatedGenes++;

    // Processing every EC number associated with the current gene.
    for (const code of codes) {

      // Splitting the EC number into its four numerical components as strings.
      const parts = code.split(".");

      // Constructing and storing the appropriate EC prefix at each classification level.
      ecLevels.forEach(({ depth }, i) => {
        byLevel[i].add(
          parts.slice(0, depth).join(".")
        );
      });
    }
  }

  // Returning the four EC profiles and the number of annotated genes.
  return {
    byLevel,
    annotatedGenes
  };
}

// Main Analysis
// Coordinating data retrieval, processing and visualization.
async function main() {

  // Loading the D3.js library.
  updateStatus("Loading D3…");
  d3 = await import("https://cdn.jsdelivr.net/npm/d3@7/+esm");

  // Retrieving pathway-gene associations.
  updateStatus("Loading pathway genes...");
  const geneRows = await runQuery("pathway_genes.sparql");
  console.log("Gene rows retrieved:", geneRows.length);

  // Retrieving gene-EC associations.
  updateStatus("Loading EC numbers...");
  const ecRows = await runQuery("enzyme_ec.sparql");
  console.log("EC rows retrieved:", ecRows.length);

  updateStatus("Calculating pathway similarities...");

  // Extracting unique pathway identifiers from the retrieved data.
  const pathwayIds = [...new Set(
    geneRows.map(row => row.wpid?.value).filter(Boolean)
  )];

  // Organizing Pathway-Gene Associations
  // Initializing a Map associating each selected pathway with an empty Set of genes.
  const genesByPathway = new Map(
    pathwayIds.map(id => [id, new Set()])
  );

  // Processing each retrieved pathway-gene association.
  for (const row of geneRows) {

    // Retrieving the gene Set associated with the pathway in the current row.
    const genes = genesByPathway.get(
      row.wpid?.value
    );

    // Adding the gene URI if the pathway and gene identifier are available.
    if (genes && row.gene?.value) {
      genes.add(row.gene.value);
    }
  }

  // Organizing Gene-EC Associations
  // Initializing a Map associating gene URIs with Sets of EC numbers.
  const ecByGene = new Map();

  // Processing each retrieved gene-EC association.
  for (const row of ecRows) {

    // Extracting the gene URI and remove surrounding whitespace from the EC number.
    const gene = row.gene?.value;
    const ec = row.ecNumber?.value?.trim();

    // Skipping rows with missing information or EC numbers that lack four numerical components.
    if (
      !gene ||
      !ec ||
      !/^\d+(?:\.\d+){3}$/.test(ec)
    ) {
      continue;
    }

    // Creating an empty EC Set if the gene has not been encountered previously.
    if (!ecByGene.has(gene)) {
      ecByGene.set(gene, new Set());
    }

    // Adding the valid EC number to the corresponding gene's Set.
    ecByGene.get(gene).add(ec);
  }

  // Comparing Pathway Pairs
  // Initialize an array to store the results of all pathway comparisons.
  const pairs = [];

  // Generating all unique pathway pairs without self-comparisons or reversed duplicates.
  for (let i = 0; i < pathwayIds.length; i++) {
    for (let j = i + 1; j < pathwayIds.length; j++) {

      // Retrieving the identifiers of the two pathways being compared.
      const idA = pathwayIds[i];
      const idB = pathwayIds[j];

      // Retrieving the Sets of genes associated with both pathways.
      const a = genesByPathway.get(idA);
      const b = genesByPathway.get(idB);

      // Checking whether both pathways have retrieved gene lists.
      const hasGeneLists =
        a.size > 0 && b.size > 0;

      // Identifing the genes unique to each pathway by removing shared genes.
      const onlyA = difference(a, b);
      const onlyB = difference(b, a);

      // Constructing EC profiles using only the nonshared genes of each pathway.
      const ecA = ecProfile(onlyA, ecByGene);
      const ecB = ecProfile(onlyB, ecByGene);

      // Checking whether both pathways have sufficient data to calculate EC similarity.
      const hasEcProfiles =
        hasGeneLists &&
        ecA.annotatedGenes > 0 &&
        ecB.annotatedGenes > 0;

      // Storing the results of the current pathway comparison.
      pairs.push({

        // Creating a label identifying the pathway pair.
        label: `${idA} vs ${idB}`,

        // Recording the number of retrieved genes in each pathway.
        geneCounts: `${a.size} / ${b.size}`,

        // Calculating the number of genes shared by both pathways.
        shared: hasGeneLists
          ? intersection(a, b).size
          : null,

        // Calculating gene overlap using the complete gene Sets.
        overlap: hasGeneLists
          ? jaccard(a, b)
          : null,

        // Recording EC annotation coverage among the nonshared genes.
        coverage: hasGeneLists
          ? `${ecA.annotatedGenes}/${onlyA.size} | ${ecB.annotatedGenes}/${onlyB.size}`
          : "Not assessable",

        // Calculating Jaccard similarity between the nonshared-gene EC profiles at all four levels.
        similarity: ecLevels.map((_, i) =>
          hasEcProfiles
            ? jaccard(ecA.byLevel[i], ecB.byLevel[i])
            : null
        )
      });
    }
  }

  // Displaying Results
  // Displaying all pathway comparisons in the web console.
  console.table(pairs);

  // Generating the four scatterplots.
  drawCharts(pairs);

  // Counting pathway pairs with assessable gene overlap and EC similarity at every level.
  const plotted = pairs.filter(
    pair =>
      pair.overlap !== null &&
      pair.similarity.every(score => score !== null)
  ).length;

  // Displaying a summary of the number of assessable pathway comparisons.
  updateStatus(
    `${plotted} pairs have gene overlap ` +
    "and EC similarity at all four levels."
  );
}

// Pearson Correlation
// Calculating Pearson r, R² and linear regression coefficients for one EC level.
function linearStatistics(points, levelIndex) {

  // Determining the number of assessable pathway pairs.
  const n = points.length;

  // Returning undefined statistics if fewer than three observations are available.
  if (n < 3) {
    return {
      n,
      r: null,
      rSquared: null,
      slope: null,
      intercept: null
    };
  }

  // Calculating the mean gene overlap.
  const xMean = d3.mean(
    points,
    pair => pair.overlap
  );

  // Calculating the mean EC similarity at the selected classification level.
  const yMean = d3.mean(
    points,
    pair => pair.similarity[levelIndex]
  );

  // Initializing the sums of squared deviations and cross-products.
  let sxx = 0;
  let syy = 0;
  let sxy = 0;

  // Calculating deviations from the means and accumulate the statistical sums.
  for (const pair of points) {
    const dx = pair.overlap - xMean;

    const dy =
      pair.similarity[levelIndex] - yMean;

    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }

  // Defining a small numerical tolerance for detecting effectively zero variation.
  const epsilon = 1e-14;

  // Returning undefined statistics if either variable has effectively no variation.
  if (sxx <= epsilon || syy <= epsilon) {
    return {
      n,
      r: null,
      rSquared: null,
      slope: null,
      intercept: null
    };
  }

  // Calculating the slope and intercept of the least-squares regression line.
  const slope = sxy / sxx;

  const intercept =
    yMean - slope * xMean;

  // Calculating Pearson's correlation coefficient and constrain it between -1 and 1.
  const r = Math.max(
    -1,
    Math.min(
      1,
      sxy / Math.sqrt(sxx * syy)
    )
  );

  // Calculating the coefficient of determination by squaring Pearson's r.
  const rSquared = r * r;

  // Returning the calculated statistical results.
  return {
    n,
    r,
    rSquared,
    slope,
    intercept
  };
}

// Preparing the Visualizations
// Preparing the common dataset and layout for all four EC classification scatterplots.
function drawCharts(pairs) {

  // Retaining only pathway pairs with assessable gene overlap and all four EC similarity scores.
  const points = pairs.filter(
    pair =>
      pair.overlap !== null &&
      pair.similarity.every(score => score !== null)
  );

  // Selecting the chart container and remove any previously generated content.
  const chart = d3.select("#trial-chart");
  chart.selectAll("*").remove();

  // Displaying an explanatory message and stop plotting if no assessable pairs exist.
  if (!points.length) {
    chart.append("p").text(
      "No pathway pairs have sufficient data to be plotted."
    );

    return;
  }

  // Creating a responsive grid layout to contain the four scatterplots.
  const grid = chart.append("div")
    .style("display", "grid")
    .style(
      "grid-template-columns",
      "repeat(auto-fit, minmax(min(100%, 340px), 1fr))"
    )
    .style("gap", "16px");

  // Generating one scatterplot for each EC classification level.
  ecLevels.forEach((level, i) => {
    drawChart(grid, points, level, i);
  });
}

// Drawing Individual Scatterplots
// Create a scatterplot displaying gene overlap and EC similarity at one classification level.
function drawChart(grid, points, level, i) {

  // Create and style a separate panel for the current scatterplot.
  const panel = grid.append("section")
    .style("background-color", "white")
    .style("border", "1px solid #dce8e2")
    .style("border-radius", "12px")
    .style("padding", "14px");

  // Display the EC classification level and its name as the panel heading.
  panel.append("h3")
    .style("margin", "0 0 8px")
    .style("color", "#163b3d")
    .text(
      `EC level ${level.depth} — ${level.name}`
    );

  // Calculate the statistical results for the current EC classification level.
  const stats = linearStatistics(points, i);

  // Create a formatting function that displays numbers with three decimal places.
  const fmt = d3.format(".3f");

  // Display Pearson r, R² and the sample size above the scatterplot.
  panel.append("p")
    .style("font-size", "0.9rem")
    .style("margin", "0 0 8px")
    .text(
      stats.r === null
        ? `Pearson r: not defined · R²: not defined · n = ${stats.n}`
        : `Pearson r = ${fmt(stats.r)} · R² = ${fmt(stats.rSquared)} · n = ${stats.n}`
    );

  // Defining Scales and SVG Container
  // Map gene overlap values from 0–20% to horizontal SVG coordinates.
  const x = d3.scaleLinear()
    .domain([0, 0.2])
    .range([65, 585]);

  // Map EC similarity values from 0–100% to reversed vertical SVG coordinates.
  const y = d3.scaleLinear()
    .domain([0, 1])
    .range([305, 25]);

  // Create a formatting function that displays decimal values as percentages.
  const percent = d3.format(".0%");

  // Create a responsive and accessibly labelled SVG container for the scatterplot.
  const svg = panel.append("svg")
    .attr("viewBox", "0 0 610 355")
    .style("width", "100%")
    .style("display", "block")
    .attr("role", "img")
    .attr(
      "aria-label",
      `Gene overlap and EC level ${level.depth} similarity with linear regression`
    );

  // Draw the horizontal axis with percentage ticks from 0% to 20%.
  svg.append("g")
    .attr("transform", "translate(0,305)")
    .call(
      d3.axisBottom(x)
        .tickValues([0, 0.05, 0.10, 0.15, 0.20])
        .tickFormat(percent)
    );

  // Draw the vertical axis with percentage ticks from 0% to 100%.
  svg.append("g")
    .attr("transform", "translate(65,0)")
    .call(
      d3.axisLeft(y)
        .ticks(5)
        .tickFormat(percent)
    );

  // Add the horizontal axis label describing gene overlap.
  svg.append("text")
    .attr("x", 325)
    .attr("y", 345)
    .attr("text-anchor", "middle")
    .attr("fill", "#273841")
    .text("Gene overlap (Jaccard index)");

  // Add the rotated vertical axis label describing EC similarity.
  svg.append("text")
    .attr(
      "transform",
      "translate(17,165) rotate(-90)"
    )
    .attr("text-anchor", "middle")
    .attr("fill", "#273841")
    .text("EC similarity (Jaccard index)");


  // Generate a unique identifier for the current chart's clipping region.
  const clipId = `trial-plot-clip-${i}`;

  // Define a rectangular clipping region matching the plotting area.
  svg.append("defs")
    .append("clipPath")
    .attr("id", clipId)
    .append("rect")
    .attr("x", 65)
    .attr("y", 25)
    .attr("width", 520)
    .attr("height", 280);

  // Create an SVG group that restricts plotted elements to the clipping region.
  const plot = svg.append("g")
    .attr(
      "clip-path",
      `url(#${clipId})`
    );

  // Draw the regression line only if valid regression statistics are available.
  if (stats.slope !== null) {

    // Determine the minimum and maximum observed gene overlap values.
    const xMin = d3.min(
      points,
      pair => pair.overlap
    );

    const xMax = d3.max(
      points,
      pair => pair.overlap
    );

    // Define a function that predicts EC similarity using the regression equation.
    const fittedY = value =>
      stats.intercept + stats.slope * value;

    // Draw an orange dashed regression line across the observed gene overlap range.
    plot.append("line")
      .attr("x1", x(xMin))
      .attr("y1", y(fittedY(xMin)))
      .attr("x2", x(xMax))
      .attr("y2", y(fittedY(xMax)))
      .attr("stroke", "#c46b43")
      .attr("stroke-width", 3)
      .attr("stroke-dasharray", "9 6");
  }

// Creating and styling the interactive tooltip.
const tooltip = panel.append("div")
  .style("position", "fixed")
  .style("background", "white")
  .style("border", "1px solid #167768")
  .style("border-radius", "8px")
  .style("padding", "10px")
  .style("font-size", "12px")
  .style("white-space", "pre-line")
  .style("pointer-events", "none")
  .style("z-index", "1000")
  .style("display", "none");

// Binding pathway-pair data to circles.
plot.selectAll("circle")
  .data(points)
  .join("circle")
  .attr("cx", pair => x(pair.overlap))
  .attr("cy", pair => y(pair.similarity[i]))
  .attr("r", 6)
  .attr("fill", "#167768")
  .attr("opacity", 0.8)
  .style("pointer-events", "all")
  .style("cursor", "pointer")

  // Displaying the tooltip when hovering over a circle.
  .on("mouseover", function(event, pair) {
    d3.select(this).attr("r", 9);

    tooltip
      .style("display", "block")
      .style("left", `${event.clientX + 12}px`)
      .style("top", `${event.clientY + 12}px`)
      .text(`${pair.label}
Genes A / B: ${pair.geneCounts}
Shared genes: ${pair.shared}
Gene overlap: ${percent(pair.overlap)}
EC similarity: ${percent(pair.similarity[i])}
EC coverage A | B: ${pair.coverage}`);
  })

  // Moving the tooltip with the mouse.
  .on("mousemove", event => {
    tooltip
      .style("left", `${event.clientX + 12}px`)
      .style("top", `${event.clientY + 12}px`);
  })

  // Hiding the tooltip when the mouse leaves the circle.
  .on("mouseout", function() {
    d3.select(this).attr("r", 6);
    tooltip.style("display", "none");
  });

  // Chart Legend
  // Creating and styling the container for the chart legend.
  const legend = panel.append("div")
    .style("display", "flex")
    .style("flex-wrap", "wrap")
    .style("align-items", "center")
    .style("gap", "8px 18px")
    .style("font-size", "0.8rem")
    .style("color", "#526269")
    .style("margin-top", "8px");

  // Creating a legend entry representing the pathway-pair markers.
  const dots = legend.append("span")
    .style("display", "inline-flex")
    .style("align-items", "center")
    .style("gap", "7px");

  // Drawing the green circular symbol used in the legend.
  dots.append("span")
    .style("display", "inline-block")
    .style("width", "10px")
    .style("height", "10px")
    .style("border-radius", "50%")
    .style("background-color", "#167768");

  // Adding the descriptive label for the green markers.
  dots.append("span")
    .text("Pathway pair (Jaccard indices)");

  // Adding a regression-line legend entry only if the regression line exists.
  if (stats.slope !== null) {

    // Creating a container for the regression-line legend entry.
    const line = legend.append("span")
      .style("display", "inline-flex")
      .style("align-items", "center")
      .style("gap", "7px");

    // Drawing an orange dashed symbol representing the regression line.
    line.append("span")
      .style("display", "inline-block")
      .style("width", "28px")
      .style(
        "border-top",
        "3px dashed #c46b43"
      );

    // Adding the descriptive label for the regression line.
    line.append("span")
      .text("Linear regression (best fit)");
  }

  // Displaying a short explanation of the interpretation of the Jaccard index.
  panel.append("p")
    .style("font-size", "0.8rem")
    .style("color", "#526269")
    .style("margin", "8px 0 0")
    .text(
      "Jaccard index: 0 = no shared entries; " +
      "1 = identical sets."
    );
}

// Executing the Analysis
// Executing the main function and handle any errors propagated during the analysis.
main().catch(error => {

  // Displaying the error message on the webpage.
  updateStatus(`Could not load the trial results. ${error.message}`);

  // Logging the error in the browser console for debugging.
  console.error(error);
});