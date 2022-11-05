# Presence-Occurrence Heatmap (POH) version 1

A simple clustered heatmap to display the overall occurence- or presence-based in a visual patterns using the R packages `gplots` and `RColorBrewer`.
Ensure that `_R-functions.R` is in the parent directory (or remove the '../' path if in the same directory) as `180522_Dendro_test.R` to implement the required functions.

This project is a really old example of one of my earlier works and it can be vastly improved on.
Several functions (such as `renameSH()`) are leftovers that should have been removed in this example repo.

The `column` variables (v1-v3) may be referred as 'isolates' in the script.
Each `row` represents a unique 'gene' with the observed counts (or a boolean presence using 0,1).

If `colorMode = TRUE` (as is in the example), the graph becomes an occurrence-based map (as shown in the Rplot.png example).
If  `colorMode = FALSE`, the graph becomes a boolean presence-based map (black(=1) if data is non-zero).

The legend is based on the row colouring which represents the intercept/number of columns (referred as 'isolates') that had a non-zero value (e.g. set 8 has been labelled  (coloured bright blue) with the maximal three isolates; sets 2, 6, and 7 have only one isolate and are coloured red; and the rest of the set combinations have two isolates present (coloured green)). 

A file `X.tsv` is output containing the ordered sets with two columns: col1 = the ordered position (from 1 to _#rows_), and col2 = the name of the set (rowname)

# Arguments
## analyze()
- The heatmap clusters (using `hclust()`) the occurrences observed
- Dataset = The data to be displayed
- DataCols = A vector containing the columns numbers to utilize for displaying the data (default: all the columns)
- Series    = The name to subset the data when wanted if applicable
- SeriesCol = the column in the provided Dataset which contains the value for the Series data
- colScale = Custom colour scale
- colorMode = If `TRUE`, display occurrence-based plot (colours the occurrences). If `FALSE`, displays a boolean presence-based plot with black (TRUE/1) or white (FALSE/0) colours.
- rowLabels = Self-explanatory
- transformData = Function name for transforming data (Default: absolute (abs)).
- showLegend = Self-explanatory
- Title = Self-explanatory
- clustAxis = which axis/data to cluster. Options: "column" (default), "row", "both", "none"
- orderPresence = If `TRUE`, reorders the rows based on the value from ReorderFunction in a descending order (largest value to the least). If `False`, plots data as read.
- ReorderFunction = What function is applied for each row (Default: sum)
- Ksplit = Obsolete, previously used to define the number of clusters wanted
- whiteBelowAbsOne = If `TRUE`, values < 1 (e.g. 0.5) would be coloured white
- ForceNumeric = Debugging option, if the data failed to be parsed as being numeric.
- DEBUG = Debugging option

# Citation
If you use this script in your project, please cite this repository and note any major modifications (if applicable) that were made.

There are several ways to cite a GitHub repository, which depends on the journals you are submitting to, but most follow a similar trend as the following:
```
Valliani, Moez. 2018. Presence-Occurrence Heatmap (POH) version 1. https://github.com/MoezV/POHeatmap. Date accessed: [date you obtained this repo]
```

It is highly recommended you use a [Referencing Software](https://en.wikipedia.org/wiki/Reference_software) if you are writing a manuscript as it makes the citation process much easier. 
