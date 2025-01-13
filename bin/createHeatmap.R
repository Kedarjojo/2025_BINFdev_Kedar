################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, metavar="path", help="Input sample file"),
  make_option(c("-g", "--geneFunctions_file"), type="character", default=NULL, metavar="path", help="Gene Functions file."),
  make_option(c("-a", "--annoData_file"), type="character", default=NULL, metavar="path", help="Annotation Data file."),
  make_option(c("-p", "--outprefix"), type="character", default='projectID', metavar="string", help="Output prefix.")
)


opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

sampleInput=opt$input_file
geneInput=opt$geneFunctions_file
annoInput=opt$annoData_file
outprefix=opt$outprefix

testing="Y"
if (testing == "Y"){
  sampleInput="sampleData.csv"
  geneInput="geneFunctions.csv"
  annoInput="annoData.csv"
  outprefix="test"
}


if (is.null(sampleInput)){
  print_help(opt_parser)
  stop("Please provide an input file.", call.=FALSE)
}

################################################
################################################
## READ IN FILES##
################################################
################################################
sampleData=read.csv(sampleInput,row.names=1)
annoData=read.csv(annoInput,row.names=1)
geneFunctions=read.csv(geneInput,row.names=1)

################################################
################################################
## Set colors##
################################################
################################################
annoColors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)

################################################
################################################
## Create a basic heatmap##
################################################
################################################

pdf(paste0("basic_heatmap_", outprefix, ".pdf"), width = 12, height = 14)
pheatmap(sampleData,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Basic Gene Expression Heatmap")
dev.off()

################################################
################################################
## Create a complex heatmap##
################################################
################################################

## Generate breaks for quantiles representing low, medium, and high categories
quantile_cutoffs <- quantile(as.matrix(sampleData), probs = c(0, 1/3, 2/3, 1))

# Open a PDF file to save the complex heatmap
output_pdf <- paste0("complex_heatmap_", outprefix, ".pdf")
pdf(output_pdf, width = 12, height = 14)

# Create the heatmap with desired settings
pheatmap(
  mat = sampleData,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D",
  annotation_col = annoData,          
  annotation_row = geneFunctions,     
  annotation_colors = annoColors,     
  show_rownames = TRUE,               
  show_colnames = TRUE,               
  annotation_names_row = FALSE,       
  annotation_names_col = FALSE,       
  breaks = quantile_cutoffs,          
  color = c("#0571b0", "#8c96c6", "#ca0020"), 
  legend_breaks = c(
    mean(quantile_cutoffs[1:2]),
    mean(quantile_cutoffs[2:3]),
    mean(quantile_cutoffs[3:4])
  ),
  legend_labels = c("Low", "Medium", "High"),
  main = "Complex Gene Expression Heatmap with Sample Annotations"
)

# Close the PDF output to save the heatmap
dev.off()
