#
#
#


library(circlize)


wpath <- "S:\\U_Proteomica\\UNIDAD\\software\\MacrosRafa\\data\\Metabolomics\\PESA_Integromics\\Data\\Analysis\\03-Network\\PESA\\ChordDiagram"
#df_name <- "qq_0.tsv"

omic2n <- list(qq=0:4, mm=0:10)



for (omic in c('qq')){
  
for (n in omic2n[[omic]]){
  
df_name <- paste0(omic,'_',n,'.tsv')

plotnamec<- paste0(strsplit(df_name, '.tsv'), "_c.png")
plotnamed<- paste0(strsplit(df_name, '.tsv'), "_d.png")
plotnamedc<- paste0(strsplit(df_name, '.tsv'), "_dc.png")

df <- read.csv(paste0(wpath, '\\', df_name), sep='\t')

#df2 <- df[(df$c==0 | df$d==0) & !(df$c==0 & df$d==0), ]
df2 <- df

dfc <- df2[,c('var1', 'var2', 'c')]
colnames(dfc) <- c('var1', 'var2', 'value')

dfd <- df2[,c('var1', 'var2', 'd')]
colnames(dfd) <- c('var1', 'var2', 'value')

elems <- dfc$var1[!duplicated(dfc$var1)]

grid.col = structure(1:length(elems), names = elems)

chordDiagram(dfc, col = ifelse(dfc$value > 0, "green", "red") ,annotationTrack = c("grid","axis"), grid.col = grid.col,
             preAllocateTracks = list(track.height = max(strwidth(elems)))
             )
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5), cex=0.5)
}, bg.border = NA) # here set bg.border to NA is important
dev.copy(jpeg, paste0(wpath, '\\',plotnamec) , width=8, height=8, units='in', res=500)
dev.off()
circos.clear()

chordDiagram(dfd, col = ifelse(dfd$value > 0, "green", "red") ,annotationTrack = c("grid","axis"), grid.col = grid.col,
             preAllocateTracks = list(track.height = max(strwidth(elems)))
)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5), cex=0.5)
}, bg.border = NA) # here set bg.border to NA is important
dev.copy(jpeg, paste0(wpath, '\\',plotnamed) , width=8, height=8, units='in', res=500)
dev.off()
circos.clear()

dfdc <- dfd
dfdc$value <- dfd$value-dfc$value
chordDiagram(dfdc, col = ifelse(dfd$value > 0, "green", "red") ,annotationTrack = c("grid","axis"), grid.col = grid.col,
             preAllocateTracks = list(track.height = max(strwidth(elems)))
)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5), cex=0.5)
}, bg.border = NA) # here set bg.border to NA is important
dev.copy(jpeg, paste0(wpath, '\\',plotnamedc) , width=8, height=8, units='in', res=500)
dev.off()
circos.clear()

}
}