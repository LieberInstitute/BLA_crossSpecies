initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "UMAP", XAxis = 1L, YAxis = 2L, 
    FacetRowByColData = "species", FacetColumnByColData = "broad_celltype", 
    ColorByColumnData = "fine_celltype", ColorByFeatureNameAssay = "logcounts", 
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample", 
    SizeByColumnData = "doubletScore", TooltipColumnData = character(0), 
    FacetRowBy = "None", FacetColumnBy = "Column Data", ColorBy = "Column data", 
    ColorByDefaultColor = "#000000", ColorByFeatureName = "LINC01409", 
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
    ColorBySampleName = "1_AAACCCACAGCGTTGC-1", ColorBySampleSource = "---", 
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
    VisualBoxOpen = FALSE, VisualChoices = c("Color", "Size"), 
    ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = .5, 
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
    CustomLabels = FALSE, CustomLabelsText = "1_AAACCCACAGCGTTGC-1", 
    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Sample", 
    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
        c(2L, 14L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
    RowSelectionSource = "---", ColumnSelectionSource = "---", 
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())

################################################################################
# Settings for Reduced dimension plot 2
################################################################################

initial[["ReducedDimensionPlot2"]] <- new("ReducedDimensionPlot", Type = "UMAP", XAxis = 1L, YAxis = 2L,
    FacetRowByColData = "species", FacetColumnByColData = "broad_celltype",
    ColorByColumnData = "fine_celltype", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample",
    SizeByColumnData = "doubletScore", TooltipColumnData = character(0),
    FacetRowBy = "None", FacetColumnBy = "Column Data", ColorBy = "Feature name",
    ColorByDefaultColor = "#000000", ColorByFeatureName = "SLC17A7",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "1_AAACCCACAGCGTTGC-1", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
    ContourColor = "#0000FF", PointSize = .5, PointAlpha = 1,
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
    CustomLabelsText = "1_AAACCCACAGCGTTGC-1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = 2L, PanelHeight = 450L, PanelWidth = 6L,
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "MIR1302-2HG", Search = "", SearchColumns = c("", 
"", "", "", "", "", ""), HiddenColumns = character(0), VersionInfo = list(
    iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version", 
    "numeric_version"))), PanelId = 1L, PanelHeight = 500L, PanelWidth = 4L, 
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
    XAxisColumnData = "fine_celltype", XAxisFeatureName = "LINC01409", 
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
    YAxisFeatureName = "SLC17A7", YAxisFeatureSource = "---", YAxisFeatureDynamicSource = FALSE, 
    FacetRowByColData = "species", FacetColumnByColData = "broad_celltype", 
    ColorByColumnData = "fine_celltype", ColorByFeatureNameAssay = "logcounts", 
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample", 
    SizeByColumnData = "doubletScore", TooltipColumnData = character(0), 
    FacetRowBy = "Column data", FacetColumnBy = "None", ColorBy = "Column data", 
    ColorByDefaultColor = "#000000", ColorByFeatureName = "LINC01409", 
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
    ColorBySampleName = "1_AAACCCACAGCGTTGC-1", ColorBySampleSource = "---", 
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, 
    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
    CustomLabelsText = "1_AAACCCACAGCGTTGC-1", FontSize = 1, 
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000", 
    VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version", 
    "numeric_version"))), PanelId = 1L, PanelHeight = 1000L, PanelWidth = 8L, 
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())
