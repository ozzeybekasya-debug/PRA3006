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

// 3. SET OPERATIONS AND JACCARD SIMILARITY
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

// 4. BUILD EC PROFILES
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

// 5. MAIN ANALYSIS
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

  // 5.1 ORGANIZE PATHWAY-GENE ASSOCIATIONS
  // Initialize a Map associating each selected pathway with an empty Set of genes.
  const genesByPathway = new Map(
    pathwayIds.map(id => [id, new Set()])
  );

  // Process each retrieved pathway-gene association.
  for (const row of geneRows) {

    // Retrieve the gene Set associated with the pathway in the current row.
    const genes = genesByPathway.get(
      row.wpid?.value
    );

    // Add the gene URI if the pathway and gene identifier are available.
    if (genes && row.gene?.value) {
      genes.add(row.gene.value);
    }
  }

  // 5.2 ORGANIZE GENE-EC ASSOCIATIONS
  // Initialize a Map associating gene URIs with Sets of EC numbers.
  const ecByGene = new Map();

  // Process each retrieved gene-EC association.
  for (const row of ecRows) {

    // Extract the gene URI and remove surrounding whitespace from the EC number.
    const gene = row.gene?.value;
    const ec = row.ecNumber?.value?.trim();

    // Skip rows with missing information or EC numbers that lack four numerical components.
    if (
      !gene ||
      !ec ||
      !/^\d+(?:\.\d+){3}$/.test(ec)
    ) {
      continue;
    }

    // Create an empty EC Set if the gene has not been encountered previously.
    if (!ecByGene.has(gene)) {
      ecByGene.set(gene, new Set());
    }

    // Add the valid EC number to the corresponding gene's Set.
    ecByGene.get(gene).add(ec);
  }

  // 5.3 COMPARE PATHWAY PAIRS
  // Initialize an array to store the results of all pathway comparisons.
  const pairs = [];

  // Generating all unique pathway pairs without self-comparisons or reversed duplicates.
  for (let i = 0; i < pathwayIds.length; i++) {
    for (let j = i + 1; j < pathwayIds.length; j++) {

      // Retrieve the identifiers of the two pathways being compared.
      const idA = pathwayIds[i];
      const idB = pathwayIds[j];

      // Retrieve the Sets of genes associated with both pathways.
      const a = genesByPathway.get(idA);
      const b = genesByPathway.get(idB);

      // Check whether both pathways have retrieved gene lists.
      const hasGeneLists =
        a.size > 0 && b.size > 0;

      // Identify the genes unique to each pathway by removing shared genes.
      const onlyA = difference(a, b);
      const onlyB = difference(b, a);

      // Construct EC profiles using only the nonshared genes of each pathway.
      const ecA = ecProfile(onlyA, ecByGene);
      const ecB = ecProfile(onlyB, ecByGene);

      // Check whether both pathways have sufficient data to calculate EC similarity.
      const hasEcProfiles =
        hasGeneLists &&
        ecA.annotatedGenes > 0 &&
        ecB.annotatedGenes > 0;

      // Store the results of the current pathway comparison.
      pairs.push({

        // Create a label identifying the pathway pair.
        label: `${idA} vs ${idB}`,

        // Record the number of retrieved genes in each pathway.
        geneCounts: `${a.size} / ${b.size}`,

        // Calculate the number of genes shared by both pathways.
        shared: hasGeneLists
          ? intersection(a, b).size
          : null,

        // Calculate gene overlap using the complete gene Sets.
        overlap: hasGeneLists
          ? jaccard(a, b)
          : null,

        // Record EC annotation coverage among the nonshared genes.
        coverage: hasGeneLists
          ? `${ecA.annotatedGenes}/${onlyA.size} | ${ecB.annotatedGenes}/${onlyB.size}`
          : "Not assessable",

        // Calculate Jaccard similarity between the nonshared-gene EC profiles at all four levels.
        similarity: ecLevels.map((_, i) =>
          hasEcProfiles
            ? jaccard(ecA.byLevel[i], ecB.byLevel[i])
            : null
        )
      });
    }
  }

  // 5.4 DISPLAY ANALYSIS RESULTS
  // Displaying all pathway comparisons in the web console.
  console.table(pairs);

  // Generating the four scatterplots.
  drawCharts(pairs);

  // Count pathway pairs with assessable gene overlap and EC similarity at every level.
  const plotted = pairs.filter(
    pair =>
      pair.overlap !== null &&
      pair.similarity.every(score => score !== null)
  ).length;

  // Display a summary of the number of assessable pathway comparisons.
  updateStatus(
    `${plotted} of ${pairs.length} pairs have gene overlap ` +
    "and EC similarity at all four levels."
  );
}

// 6. PEARSON CORRELATION AND LINEAR REGRESSION
// Calculate Pearson r, R² and linear regression coefficients for one EC level.
function linearStatistics(points, levelIndex) {

  // Determine the number of assessable pathway pairs.
  const n = points.length;

  // Return undefined statistics if fewer than three observations are available.
  if (n < 3) {
    return {
      n,
      r: null,
      rSquared: null,
      slope: null,
      intercept: null
    };
  }

  // Calculate the mean gene overlap.
  const xMean = d3.mean(
    points,
    pair => pair.overlap
  );

  // Calculate the mean EC similarity at the selected classification level.
  const yMean = d3.mean(
    points,
    pair => pair.similarity[levelIndex]
  );

  // Initialize the sums of squared deviations and cross-products.
  let sxx = 0;
  let syy = 0;
  let sxy = 0;

  // Calculate deviations from the means and accumulate the statistical sums.
  for (const pair of points) {
    const dx = pair.overlap - xMean;

    const dy =
      pair.similarity[levelIndex] - yMean;

    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }

  // Define a small numerical tolerance for detecting effectively zero variation.
  const epsilon = 1e-14;

  // Return undefined statistics if either variable has effectively no variation.
  if (sxx <= epsilon || syy <= epsilon) {
    return {
      n,
      r: null,
      rSquared: null,
      slope: null,
      intercept: null
    };
  }

  // Calculate the slope and intercept of the least-squares regression line.
  const slope = sxy / sxx;

  const intercept =
    yMean - slope * xMean;

  // Calculate Pearson's correlation coefficient and constrain it between -1 and 1.
  const r = Math.max(
    -1,
    Math.min(
      1,
      sxy / Math.sqrt(sxx * syy)
    )
  );

  // Calculate the coefficient of determination by squaring Pearson's r.
  const rSquared = r * r;

  // Return the calculated statistical results.
  return {
    n,
    r,
    rSquared,
    slope,
    intercept
  };
}

// 7. PREPARE THE FOUR SCATTERPLOTS
// Prepare the common dataset and layout for all four EC classification scatterplots.
function drawCharts(pairs) {

  // Retain only pathway pairs with assessable gene overlap and all four EC similarity scores.
  const points = pairs.filter(
    pair =>
      pair.overlap !== null &&
      pair.similarity.every(score => score !== null)
  );

  // Select the chart container and remove any previously generated content.
  const chart = d3.select("#trial-chart");
  chart.selectAll("*").remove();

  // Display an explanatory message and stop plotting if no assessable pairs exist.
  if (!points.length) {
    chart.append("p").text(
      "No pathway pairs have sufficient data to be plotted."
    );

    return;
  }

  // Create a responsive grid layout to contain the four scatterplots.
  const grid = chart.append("div")
    .style("display", "grid")
    .style(
      "grid-template-columns",
      "repeat(auto-fit, minmax(min(100%, 340px), 1fr))"
    )
    .style("gap", "16px");

  // Generate one scatterplot for each EC classification level.
  ecLevels.forEach((level, i) => {
    drawChart(grid, points, level, i);
  });
}

// 8. DRAW EACH INDIVIDUAL SCATTERPLOT
// Create a scatterplot displaying gene overlap and EC similarity at one classification level.
function drawChart(grid, points, level, i) {

  // 8.1 CREATE THE CHART PANEL
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

  // 8.2 DISPLAY STATISTICAL RESULTS
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

  // 8.3 DEFINE SCALES AND SVG CONTAINER
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

  // 8.4 CREATE AXES AND AXIS LABELS
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


  // 8.5 DEFINE THE VISIBLE PLOTTING AREA
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


  // 8.6 DRAW THE LINEAR REGRESSION LINE
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


  // 8.7 DRAW THE PATHWAY PAIRS
  // Bind pathway-pair data to green circles positioned by gene overlap and EC similarity.
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

    // Adding a tooltip displaying detailed pathway comparison results.
    .append("title")
    .text(pair => `${pair.label}
Genes A / B: ${pair.geneCounts}
Shared genes: ${pair.shared}
Gene overlap: ${percent(pair.overlap)}
EC level ${level.depth} similarity: ${percent(pair.similarity[i])}
EC coverage A | B: ${pair.coverage}`);

  // 8.8 CREATE THE CHART LEGEND
  // Create and style the container for the chart legend.
  const legend = panel.append("div")
    .style("display", "flex")
    .style("flex-wrap", "wrap")
    .style("align-items", "center")
    .style("gap", "8px 18px")
    .style("font-size", "0.8rem")
    .style("color", "#526269")
    .style("margin-top", "8px");

  // Create a legend entry representing the pathway-pair markers.
  const dots = legend.append("span")
    .style("display", "inline-flex")
    .style("align-items", "center")
    .style("gap", "7px");

  // Draw the green circular symbol used in the legend.
  dots.append("span")
    .style("display", "inline-block")
    .style("width", "10px")
    .style("height", "10px")
    .style("border-radius", "50%")
    .style("background-color", "#167768");

  // Add the descriptive label for the green markers.
  dots.append("span")
    .text("Pathway pair (Jaccard indices)");

  // Add a regression-line legend entry only if the regression line exists.
  if (stats.slope !== null) {

    // Create a container for the regression-line legend entry.
    const line = legend.append("span")
      .style("display", "inline-flex")
      .style("align-items", "center")
      .style("gap", "7px");

    // Draw an orange dashed symbol representing the regression line.
    line.append("span")
      .style("display", "inline-block")
      .style("width", "28px")
      .style(
        "border-top",
        "3px dashed #c46b43"
      );

    // Add the descriptive label for the regression line.
    line.append("span")
      .text("Linear regression (best fit)");
  }

  // Display a short explanation of the interpretation of the Jaccard index.
  panel.append("p")
    .style("font-size", "0.8rem")
    .style("color", "#526269")
    .style("margin", "8px 0 0")
    .text(
      "Jaccard index: 0 = no shared entries; " +
      "1 = identical sets."
    );
}

// 10. EXECUTE THE ANALYSIS
// Execute the main function and handle any errors propagated during the analysis.
main().catch(error => {

  // Display the error message on the webpage.
  updateStatus(`Could not load the trial results. ${error.message}`);

  // Log the error in the browser console for debugging.
  console.error(error);
});
