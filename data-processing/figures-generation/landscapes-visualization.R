###
## An example of landscapes visualization
## Argument x of heatmap.3 should be a landscape matrix to be visualized
###

psz(10)
#pdf(file="heatmap_landscape_manhattan_wardD2_final.pdf", width = 10, height=10 )
par(cex.axis = 1.5, mar=c(0,0,0,0))
#png(file="heatmap_landscape_manhattan_wardD2_nodend.png", 1000, 1000)
heatmap.3(
        x = as.matrix(kidney.data$landscape$landscapes[non.constant>0.0001,]), 
        KeyValueName="Activity", 
        labCol = F, labRow = F, 
        margins = c(1.75,1.75), 
        lwid = c(1.3,0.15,5),
        lmat=rbind(c(6,0,5), c(0,0,2), c(4,1,3)), 
        xlab="Metabolic Landscape", ylab="Reaction Activity",
        lhei=c(0.75,0.25,3.5),  
        ColSideColors = colSide, 
        ColSideColorsSize = 3, col=gray.colors(10, 0, 1),
        
        RowSideColorsSize = 1.5, 
        RowSideColors = rowSide.land,# matrix(nrow = 1, 
        #key=TRUE, symkey=FALSE, symbreaks=FALSE, keysize=1, 
        hclustfun=myclust, distfun=mydist, 
        #density.info="none", trace="none", 
        )
        legend(x = 0, y= 0.85,x.intersp = 0.3, 
          legend=c("","Tissue Type (TT)",
                   "Primary Tumor","Solid Normal",
                   "","Morph. Type (MT)",
                   "Clear Cell RCC","Papillary RCC","Chromophobe RCC", "Unclassified",
                   "","Gene Rules (R)",
                   "SLCO1B1", "SLCO1A2", "SLC7A5", "SLC7A9", "other"),
          fill=c("white","white",
                 "red2","green2",
                 "white","white",
                 "yellow2","orange2","blue2", "lightblue",
                 "white","white",
                 "chocolate1", "coral1","goldenrod2","darksalmon", "grey50"), 
          border=FALSE, bty="n",y.intersp = 1, cex=1)
#dev.off()
