library(circlize)

wpath <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\03-Network\\PESA\\ChordDiagram"

df_name <- "meta-network_WGCNA.tsv"
plotname <- "meta-network_WGCNA.png"

df <- read.csv(paste0(wpath, '\\', df_name), sep='\t')

#elems <- union(df$var1, df$var2)
elems <- df$var1[!duplicated(df$var1)]

grid.col1 = structure(1:length(elems), names = elems)

elems <- df$var2[!duplicated(df$var2)]

grid.col2 = structure(rep('grey', length(elems)), names = elems)

grid.col = c(grid.col1, grid.col2)


chordDiagram(df ,annotationTrack = c("grid","axis"), grid.col = grid.col,
             preAllocateTracks = list(track.height = max(strwidth(elems)))
)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(-1, 0.5), cex=1)
}, bg.border = NA) # here set bg.border to NA is important
dev.copy(jpeg, paste0(wpath, '\\',plotname) , width=8, height=8, units='in', res=500)
dev.off()
circos.clear()
