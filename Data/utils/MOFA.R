library(MOFA2)
library(ggplot2)

xq_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Proteomics\\AWHS\\WorkingFiles\\Xq_minus_X_norm.tsv"
xm_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Metabolomics\\AWHS\\WorkingFiles\\Xm_norm_MS2.tsv"

mdata_path <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Metadata\\AWHS\\WorkingFiles\\main_metadata.tsv"

xq <- read.table(xq_path, sep="\t", header=T, row.names=1)
xm <- read.table(xm_path, sep="\t", header=T, row.names=1)


seqn <- intersect(
  rownames(xq),
  rownames(xm)
)

mdata <- read.table(mdata_path, sep="\t", header=T, row.names=1)[seqn,]

data <- list(
  q=t(as.matrix(xq[seqn,])),
  m=t(as.matrix(xm[seqn,]))
)

MOFAobject <- create_mofa(data)

plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 20
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path('S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\04-FactorAnalysis\\AWHS',"model.hdf5")
model <- run_mofa(MOFAobject, outfile)


#
# Exploration
#

mdata['sample'] <- rownames(mdata)
samples_metadata(model) <- mdata

head(model@cache$variance_explained$r2_total[[1]])
plot_variance_explained(model, x="view", y="factor")
plot_variance_explained(model, x="group", y="factor", plot_total = T)

plot_factor(model, 
  factor = 1:5,
  color_by = "Caso.control",
  #shape_by = "condition"
)


p <- plot_factor(model, 
                 factors = 1:5,
                 color_by = "Caso.control",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
#p <- p + 
#  scale_color_manual(values=c('0'="black", '1'="red")) +
#  scale_fill_manual(values=c('0'="black", '1'="red"))

print(p)


plot_factors(model,
             factors = 1:5,
             color_by = "Caso.control"
)


plot_data_heatmap(model,
                  view = "m",         # view of interest
                  factor = 4,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE
)

dimr <- run_umap(model)
plot_dimred(model,
            method = "UMAP",  # method can be either "TSNE" or "UMAP"
            color_by = "Caso.control"
)


pdf(file=file.path('S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\04-FactorAnalysis\\AWHS',"plots.pdf"))
dev.off()
