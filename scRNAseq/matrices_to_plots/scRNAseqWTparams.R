# We define the color palette used for featurePlots 
pal <- c("#FEB24C",	"#FD8D3C",	"#FC4E2A",	"#E31A1C",	"#BD0026",	"#800026")
Col.featurePlot <- alpha(c("#D3D3D3", pal),0.9)

# And the color palette for UMAP:

# The time:
my.time.colors <- grey(level = c(0.9,0.85,0.7,0.55,0.4,0.25,0.1))
names(my.time.colors) <- c("0h", "48h", "72h", "96h", "120h", "144h", "168h")

# And the clusters:
# Group new cluster names by categories so they will have colors along a colormap:
# list.Fate.level is a list with the categories:
# Order matter! It will be used in all plots!
list.Fate.level <- list("Pluripotent" = c("ESCs", "Epiblast", "Pluripotent" ),
                        "PrimStreak" = c("Prim. Streak", "Ant. Prim. Streak", "Caudal Mes."),
                        "Neuronal" = c("Early NMP","Late NMP", "Neural Tube", "Neuron Progen. 1", "Neuron Progen. 2","Neuron Precursor"),
                        "Mesoderm" = c("Post. PSM", "Mixed Mes.", "Ant. PSM", "Som. Mes.", "Dermomyotome", "Sclerotome"),
                        "Endoderm" = c("Early Endoderm", "Late Endoderm"),
                        "Other" = c("Endothelial", "Surface Ecto.", "Unknown"))
# list.color is a list with the same categories and colors that describe the heatmap. There must be at least 2 colors:
list.color  <- list("Pluripotent" = c('#000000' , '#BFBFBF'),
                    "PrimStreak" =c('#400000','#FFF2F2'),
                    "Neuronal" = c('#DAE3F3', '#002060'),
                    "Mesoderm" = c('#FBE5D6','#5F2C09'),
                    "Endoderm" = c('#5FE756', '#70AD47'),
                    "Other" = c('#FBBEDE', '#7030A0'))
# Check all fates have a category name
if ("" %in% names(list.Fate.level)) {
  stop("Some fates level has no name!\n")
}
# Check all category names in list.Fate.level are also in list.color
if (any(! names(list.Fate.level) %in% names(list.color))) {
  stop("The following Fates are not in colors: ", paste(names(list.Fate.level)[!names(list.Fate.level) %in% names(list.color)], collapse = ", "), "\n")
}
# Create the colormaps
my.fate.colors <- unlist(lapply(names(list.Fate.level), function(fate) {
  colorRampPalette(list.color[[fate]])(length(list.Fate.level[[fate]]))
}))
# Give them as names the new cluster names
names(my.fate.colors) <- unlist(list.Fate.level)