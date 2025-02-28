library(Seurat) # Version 4 as in MS
library(dplyr)
library(gdata)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(SummarizedExperiment)
library(BiocParallel)
library(xlsx)
library(sctransform)
library(scRepertoire)
library(data.table)
library(magrittr)
library(monocle3)
library(tidyr)
library(patchwork)
library(magrittr)
library(SeuratWrappers)
library(viridis)
library(pheatmap)
set.seed(1000)



#Read filtered 10x count data for each sample  
nraw_counts <- Read10X(data.dir = "...")
n = "Ad_Non"

summary(colSums(nraw_counts))

# check out the first six genes and cells
nraw_counts[1:10, 1:20]


# check how many genes have at least one transcript in each cell
at_least_one <- apply(nraw_counts, 2, function(x) sum(x > 0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(nraw_counts),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")

tmp <- apply(nraw_counts, 1, function(x) sum(x>0))
table(tmp>=5) # Number of genes in 5 or more cells



## Create seurat object
s1_SO<- CreateSeuratObject(counts = nraw_counts, min.cells = 10, min.features = 200, project = paste0(n,"_10x_scRNAseq"))
s1_SO[["percent.mt"]] <- PercentageFeatureSet(s1_SO, pattern = "mt-")
VlnPlot(s1_SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

quantile(apply(GetAssayData(s1_SO,slot = "counts"), 2, function(x) sum(x > 0)), c(.25, .50,  .75, .90, .99))
quantile(apply(GetAssayData(s1_SO,slot = "counts"), 2, function(x) sum(x)), c(.25, .50,  .75, .90, .99))

##wb131 cutoff 6154
#total cells 1220

##wb132 cutoff 4659
#total cells 1266

##wb133 cutoff 5286
#total cells 1559

#wb134 cutoff 3542
#total cells 1978


s1_SO <- subset(s1_SO, subset =  nFeature_RNA < 3542 &  nCount_RNA > 1000 & percent.mt < 5)
ncol(GetAssayData(s1_SO,slot = "counts"))
summary(apply(GetAssayData(s1_SO,slot = "counts"), 2, function(x) sum(x>0)))



n1_pbmc <- SCTransform(s1_SO, vst.flavor = "v2", vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
n2_pbmc <- SCTransform(s2_SO, vst.flavor = "v2", vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
n3_pbmc <- SCTransform(s3_SO, vst.flavor = "v2", vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
n4_pbmc <- SCTransform(s4_SO, vst.flavor = "v2", vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)

#SavePBMC

#################################
######### Subclustering #########
#################################

s_SO <- s4_SO
n=4
n_cd4 <- subset(s_SO, Cd4 > 0)#
n_cd4 <- subset(n_cd4, Cd8a == 0)

n_cd4c <- Cells(n_cd4)
t <- subset(s_SO, cells = n_cd4c)#
ncol(GetAssayData(t,slot = "counts"))
n_cd4r <- GetAssayData(t, slot = "counts")


s_SO<- CreateSeuratObject(counts = n_cd4r, project = paste0("S",n,"_Cd4"))
s_SO[["percent.mt"]] <- PercentageFeatureSet(s_SO, pattern = "mt-")#
n4cd4_pbmc <- SCTransform(s_SO, vst.flavor = "v2", vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)#

n1cd8_pbmc$stim <- "Neo_Inf"
n2cd8_pbmc$stim <- "Neo_Non"
n3cd8_pbmc$stim <- "Ad_Inf"
n4cd8_pbmc$stim <- "Ad_Non"

##
n1_pbmc$stim <- "Neo_Inf"
n2_pbmc$stim <- "Neo_Non"
n3_pbmc$stim <- "Ad_Inf"
n4_pbmc$stim <- "Ad_Non"


##Integration_standard
immune.anchors <- FindIntegrationAnchors(object.list = list(n1_pbmc,n3_pbmc), dims = 1:30)
pbmc.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
DefaultAssay(pbmc.combined) <- "integrated"
pbmc <- ScaleData(pbmc.combined, verbose = FALSE)

###################
##Integration_SCT##
###################

options(future.globals.maxSize = 4000 * 1024^2)

n.features <- SelectIntegrationFeatures(object.list = list(n1cd8_pbmc,n2cd8_pbmc,n3cd8_pbmc,n4cd8_pbmc), nfeatures = 3000)
n.list <- PrepSCTIntegration(object.list = list(n1cd8_pbmc,n2cd8_pbmc,n3cd8_pbmc,n4cd8_pbmc), anchor.features = n.features, 
                             verbose = FALSE)
n.anchors <- FindIntegrationAnchors(object.list = n.list, normalization.method = "SCT", 
                                    anchor.features = n.features, verbose = FALSE)
n.integrated <- IntegrateData(anchorset = n.anchors, normalization.method = "SCT", 
                              verbose = FALSE)
marrow <- RunPCA(n.integrated, verbose = FALSE)

"clusternormally"

DefaultAssay(pbmc2) <- "RNA"
pbmc.integrated <- NormalizeData(pbmc2, verbose = FALSE)


##########################################


marrow <- RunPCA(pbmc,verbose = FALSE)

DimPlot(marrow)


pbmc2 <- RunUMAP(marrow, dims = 1:30, verbose = FALSE)

pbmc2 <- FindNeighbors(pbmc2, dims = 1:30, verbose = FALSE)
pbmc2 <- FindClusters(pbmc2, verbose = FALSE)


## MS
tiff("All_UMAP_ms2.tiff",units="cm", width=10, height=8, res=300)
DimPlot(pbmc.integrated, group.by = "stim", pt.size = 0.5) & scale_color_manual(values = c('#e95462',"black")) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),aspect.ratio = 1) & xlim(-6.5,6.5) & ylim(-6.5,6.5) &
  guides(color = guide_legend(override.aes = list(size=1), ncol=2)) & labs(title = "Inf vs Non-Inf")
dev.off()


#MS
tiff("All_UMAP.tiff",units="cm", width=10, height=8, res=300)
DimPlot(pbmc.integrated, pt.size = 0.5) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),aspect.ratio = 1) & xlim(-6,8) & ylim(-6,8) &
  guides(color = guide_legend(override.aes = list(size=2,shape =18), keyheight = 0.2,keywidth = 0.5, ncol=2)) & labs(title = "All samples")
dev.off()
#DimPlot(pbmc2, split.by = "stim")
#DimPlot(pbmc2, group.by = "stim")
#DimPlot(pbmc2, label = TRUE, group.by ="orig.ident") 

##MS
#cd4  c("#56B4E9",'#de3230',"#ff7321","black",'#f6d746')
tiff("Cd8a/Cd8a_combined-3.tiff",units="cm", width=10, height=8, res=300)
DimPlot(cd8_custom, pt.size = 0.2, cols = c("#56B4E9",'#de3230','#f6d746',"grey")) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),aspect.ratio = 1) & xlim(-6,6) & ylim(-6,6) &
  guides(color = guide_legend(override.aes = list(size=1), ncol=2)) & labs(title = "Cd8a") 
dev.off()

#MS
tiff("Cd4/Cd4_genes-1.tiff",units="in", width=10, height=10, res=420)
FeaturePlot(pbmc.integrated,pt.size = 0.13,features = c("Foxp3","Il2ra","Cxcr6","Ly6c2","Cd44","Tcf7","Ccr7","Sell","Mki67","Pclaf","Cxcr5","Izumo1r"), ncol = 3) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.line =element_line(size = 0.5)) & xlim(-6,8) & ylim(-6,8) &
  theme(legend.key.size = unit(0.2, "cm"), aspect.ratio = 1)
dev.off()

#Cytotoxicity

genes <- c("Gzmm", "Fasl", "Gzmc", "Gzmb", "Gzma", "Gzmk", "Prf1", "Ifng", "Tnf", "Tnfsf10")

tiff("Cd8a/Cd8_module-score_regulation.tiff",units="in", width=8, height=5, res=220)
FeaturePlot(pbmc1,pt.size = 0.2,  features = "regulation1", split.by = "stim")  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1.5,1.5)) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.line =element_line(size = 0.5)) & xlim(-6,8) & ylim(-6,8) &
  theme(legend.position = c(0.9,0.87),legend.key.size = unit(0.20, "cm"), aspect.ratio = 1) 
dev.off()

##
tiff("Cd8a/Cd8_module-score_regulation_b_violin.tiff",units="in", width=6, height=4, res=120)
VlnPlot(subset(pbmc1,subset = stim %in% c("Ad_Inf","Neo_Inf")),pt.size = 0.8,features = "regulation_b1",split.by = "stim",cols = c("#56B4E9",'#de3230'))
dev.off()
##

## "Naive"   "Th1"     "Tcmp"    "Treg"    "Cycling"
tiff("Cd4/Cd4_module-score_exhaustion_Cycling.tiff",units="in", width=5, height=4, res=120)
VlnPlot(subset(pbmc1,subset = seurat_clusters2 == "Cycling" & stim %in% c("Ad_Inf","Neo_Inf")),pt.size = 0.8,split.by = "stim",features = "exhaustion1",cols = c("#56B4E9",'#de3230')) + stat_compare_means(method = "wilcox")
dev.off()

tiff("Cd4/Cd4_module-score_exhaustion.tiff",units="in", width=8, height=5, res=220)
FeaturePlot(pbmc1,pt.size = 0.2,  features = "exhaustion1", split.by = "stim")  & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-1.5,1.5)) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.line =element_line(size = 0.5)) & xlim(-6,8) & ylim(-6,8) &
  theme(legend.position = c(0.9,0.87),legend.key.size = unit(0.20, "cm"), aspect.ratio = 1) 
dev.off()



#MS  
tiff("PCA_all_samples2.tiff", units="in", width=5.2, height=4.2, res=150)
ggplot(dataGG, aes(x=PC1,y=PC2, label = rownames(dataGG), color=condition, shape=sample2)) +  
  geom_point(size=3.1, stroke = 2) +   geom_text_repel(size=3.15,stat = "identity") +
  theme(axis.text.x  = element_text(size = 4.5,margin = margin(r = 0))) +
  theme(axis.text.y = element_text(size = 4.5,margin = margin(r = 0))) +
  xlab(paste0("PC1, VarExp:", round(percentVar[1],4))) + 
  ylab(paste0("PC2, VarExp:", round(percentVar[2],4))) + 
  theme(axis.title.y = element_text(size = 10.5))+
  theme(axis.title.x = element_text(size = 10.5))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2), aspect.ratio = 1) + 
  scale_color_manual(values=c("#e95462","black"))  + scale_shape_manual(values = c(19,1)) +
  theme(legend.position="none") + ylim(-40,40) + xlim(-40,40)
dev.off()


#c('#e95462','#22a884','#f6d746','#331067')
#c('grey','grey','grey','#331067')  
#  scale_shape_manual(values = c(20,17,1,8))

## Dotplot MS
tiff(paste0("./Cd8a/Non_Inf_dotplot.tiff"),units="in", width=7, height=5, res=130)
ggplot(meta_summary2, aes(x=Gene, y= interaction)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21,stroke = 0.1) +
  scale_size("% detected", range = c(1,11)) +
  scale_fill_gradientn(colours = viridis::plasma(100, direction = 1),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression") +
  ylab("Cluster") + xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size=9, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title = element_text(size=7), aspect.ratio = 0.8) + scale_y_discrete(limits=rev) 
dev.off()

## Density MS
tiff(paste0("./Cd4/Neo_Inf-1.tiff"),units="in", width=10, height=8, res=130)
ggplot(data2[data2$cluster=="Neo_Inf",]) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color =  cluster), size=2.3, show.legend = TRUE) +
  stat_density_2d(aes(x=UMAP_1, y=UMAP_2,alpha = (..level..)), colour ='#e95462', bins =5, size = 1.5) +  
  scale_color_manual(values=c('#e95462')) + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),aspect.ratio = 1) +
  
  ylim(-6,8) + xlim(-6,8) + 
  theme(text = element_text(size = 8),legend.key.size = unit(0.06, 'cm')) +  scale_alpha_continuous(range = c(0, 1))
dev.off()



## Dotplot MS
genes <- c("Sell","Ccr7","Tcf7","Mki67","Pclaf","Pdcd1","Cd44","Ccl5","Cxcr6","Foxp3","Lag3","Tbx21","Cxcr5",  
           "Izumo1r","Il10","Il21","Ly6c2","Ly6a")
genes <- Hmisc::Cs(Ccr7, Tcf7, Il7r, Sell, Ly6c2, Cd44, Gzmb, Gzmk, Ccl5, Ccl3, Xcl1, Pclaf, Mki67, Top2a, Entpd1, Havcr2, Tigit, Lag3, Pdcd1, Tox)
genes <- c("Ccr7", "Tcf7", "Sell", "Mki67", "Pclaf", "Xcl1", "Gzmb", "Gzmk", "Ccl5", "Cxcr6", "S1pr5", "Klrg1", "Cx3cr1", "Itga1")
genes <- c("Sell", "Ms4a4c", "Satb1", "Ctsd", "Itgb1", "S100a10", "Lgals3", "Gzma", "Ifitm1", "Gzmk", "Itga4", "Ccl5")
genes <- Hmisc::Cs(Ctla2a, Il2rb, Klra7, Nsg2, Txk, Eomes, Gzma, Gzmk, Lgals3, Itga4, Mki67, Isg15, Ifit3, Ly6a, Gzmb, Zbp1, Cxcl10, Klf2, Klf3, Lef1, Tsc22d3, 
                   Il7r, S100a6, Ifitm1, Cxcr6, Id2,Cxcr3, Ly6c2, Ccr2, Nkg7, Ccr7, Dapl1, Sell, Ccl5, S1pr5, Cx3cr1, Klrg1, Zeb2, Itgb2, Lgals1, Tbx21, Tcf7, Ly6e, Lef1, Plac8)

exp_mat <- as.matrix(cd8_custom2@assays$RNA[genes,])
meta <- cd8_custom2@meta.data %>% 
  select(seurat_clusters)
meta <- bind_cols(meta,as.data.frame(cd8_custom2$stim), as.data.frame(t(exp_mat)))
colnames(meta)[2] <- "stim"

meta <- pivot_longer(meta, c(-seurat_clusters,-stim), names_to="Gene", values_to="Expression")
meta_summary <- meta %>%
  group_by(seurat_clusters,stim, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)

meta_summary2 <- meta_summary[order(meta_summary$stim,meta_summary$seurat_clusters,meta_summary$Avg),]

meta_summary2$Gene <- factor(meta_summary2$Gene,levels =  genes)

#c("Naive","Cycling","TCM","Effector")
#Naive Cycling Tcmp Th1 Treg  
#"Ad_Non","Neo_Non", "Ad_Inf",  "Neo_Inf" 

meta_summary2$stim <- factor(meta_summary2$stim, levels = c( "Ad_Non","Neo_Non", "Ad_Inf",  "Neo_Inf" ))

meta_summary2$seurat_clusters <- factor(meta_summary2$seurat_clusters, levels = c("Naive","Cycling","TCM","Effector"))

meta_summary2 <- data.frame(with(meta_summary2, meta_summary2[order(seurat_clusters,stim),]))
meta_summary2$interaction <- interaction(meta_summary2$stim,meta_summary2$seurat_clusters)

tiff("./Cd4/Cd4_dotplot.tiff", units="cm", width=36, height=24, res=140)

ggplot(meta_summary2, aes(x=Gene, y= interaction)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21,stroke = 0.1) +
  scale_size("% detected", range = c(1,11)) +
  scale_fill_gradientn(colours = viridis::plasma(100, direction = 1),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression") +
  ylab("Cluster") + xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size=9, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title = element_text(size=7), aspect.ratio = 1) + scale_y_discrete(limits=rev) 

dev.off()

#MS
tiff(paste0("./Cd8a/Inf_Effector_dotplot.tiff"),units="in", width=7, height=3, res=130)
ggplot(meta_summary2[grep("Inf.Eff",meta_summary2$interaction),], aes(x=Gene, y= interaction)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21,stroke = 0.1) +
  scale_size("% detected", range = c(1,11)) +
  scale_fill_gradientn(colours = viridis::plasma(100, direction = 1),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression") +
  ylab("Cluster") + xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size=9, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title = element_text(size=7), aspect.ratio = 0.4) + scale_y_discrete(limits=rev) 
dev.off()




pheatmap(test1[grep("1_3",colnames(test1))], cluster_rows=TRUE, show_rownames=TRUE,cluster_cols=FALSE,scale = "none", 
         show_colnames = FALSE, col=viridis::plasma(25, direction = 1),annotation_col=df[grep("1_3",rownames(df[1.])),][1.])

##PCA
p1 <-   DimPlot(cd8_combined, group.by = "stim", pt.size = 1.7, cols = c('#331067','grey','grey','grey'), shape.by = "stim", order =  FALSE) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank()) & scale_shape_manual(values = c(21,21,21,21))
p2 <-   DimPlot(cd8_combined, group.by = "stim", pt.size = 1.7, cols = c('grey','#331067','grey','grey'), shape.by = "stim", order =  FALSE) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank())& scale_shape_manual(values = c(21,21,21,21))
p3 <-   DimPlot(cd8_combined, group.by = "stim", pt.size = 1.7, cols = c('grey','grey','#331067','grey'), shape.by = "stim", order =  FALSE) & 
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank())& scale_shape_manual(values = c(21,21,21,21))
p4 <-   DimPlot(cd8_combined, group.by = "stim", pt.size = 1.7, cols = c('grey','grey','grey','#331067'), shape.by = "stim", order =  FALSE)& 
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank())& scale_shape_manual(values = c(21,21,21,21))
wrap_plots(p1, p2,p3,p4)  



tiff("PCA_pseudo_bulk.tiff", units="in", width=5.2, height=4.2, res=300)
ggplot(dataGG, aes(x=PC1,y=PC2, label = rownames(dataGG), color=sample)) + geom_point(size=3) +xlab(paste0("PC1, VarExp:", round(percentVar[1],4))) + 
  ylab(paste0("PC2, VarExp:", round(percentVar[2],4)))
dev.off()


Idents(pbmc2) <- "stim"

avg.t.cells <- as.data.frame(AggregateExpression(pbmc2, slot = 'counts')$RNA)

avg.t.cells <- as.data.frame(log1p(AverageExpression(pbmc.integrated, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

tiff("GFP.tiff", units="cm", width=19, height=8, res=200)
FeaturePlot(pbmc.integrated, features = "mt-GFP", split.by = "stim",pt.size =0.4) & theme(text = element_text(size = 5),
                                                                                          axis.text.y = element_blank(),
                                                                                          axis.text.x = element_blank(),
                                                                                          axis.ticks.x=element_blank(),
                                                                                          axis.ticks.y = element_blank())
#ggplot(vsn, aes(Neo_Inf, Neo_Non)) + geom_point() +geom_cor(method = "kendall", ypos = 1e5)
dev.off()

##################### DE scRNA ##########################

DefaultAssay(pbmc2) <- "RNA"
pbmc.integrated <- NormalizeData(pbmc2, verbose = FALSE)

##subset
cluster2 <- Idents(pbmc.integrated)

Idents(pbmc.integrated_) <- "orig.ident"
Idents(pbmc.integrated_) <- cluster2

pbmc.integrated <- subset(x = pbmc.integrated_, idents = "4")

Idents(pbmc.integrated) <- rep("CD8",4)
pbmc.integrated$celltype.stim <- paste(Idents(pbmc.integrated), pbmc.integrated$stim, sep = "_")
pbmc.integrated$celltype <- Idents(pbmc.integrated)
Idents(pbmc.integrated) <- "celltype.stim"

DE_inf <- FindMarkers(pbmc.integrated, ident.1 = "CD8_Neo_Inf", ident.2 = "CD8_Ad_Inf", verbose = FALSE, min.pct = 0,logfc.threshold = 0)

DE_inf <- DE_inf[order(DE_inf$avg_log2FC,decreasing = TRUE),]

DE_inf %>% slice_max(avg_log2FC,n=20) -> DE_inf1
DE_inf %>% slice_min(avg_log2FC,n=21) -> DE_inf2

DE_inf_ <- bind_rows(DE_inf1,DE_inf2)
DE_inf_ <- DE_inf_[order(DE_inf_$avg_log2FC,decreasing = TRUE),][-1,]

write.xlsx(DE_inf[order(DE_inf$avg_log2FC,decreasing = TRUE),],"Cd8a/Ad_Neo-Inf_All-DE_cluster4.xlsx")

for (i in 1:length(rownames(DE_inf_))){
  png(paste0("Cd8a/",rownames(DE_inf_)[i],"_Cd8_Ad_Non-Inf_cluster4.png"),units="in", width=6, height=5, res=120)
  #p1 <- FeaturePlot(pbmc.integrated, features = head(rownames(DE_inf[order(DE_inf$avg_log2FC,decreasing = TRUE),]))[i], split.by = "stim",pt.size = 1.5)
  #p2 <- VlnPlot(pbmc.integrated, features = head(rownames(DE_inf[order(DE_inf$avg_log2FC,decreasing = TRUE),]))[i])
  print(VlnPlot(pbmc.integrated, features = rownames(DE_inf_)[i],pt.size = 1.5))
  #print(plot_grid(p1,p2, ncol = 1))
  dev.off()  
}


#######################################

md <- pbmc2@meta.data %>% as.data.table
md[, .N, by = c("stim", "seurat_clusters")] %>% dcast(., stim ~ seurat_clusters, value.var = "N")

prop.table(table(Idents(pbmc)))


DefaultAssay(pbmc.integrated) <- "integrated"


c.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
head(c.markers)

#MS heatmap
df <- data.frame()
for (i in 1:length(x)){
  DefaultAssay(pbmc.integrated) <- "RNA"
  c.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = x[i], grouping.var = "stim", verbose = FALSE)
  colnames(c.markers) <- gsub("__","",colnames(c.markers))
  c.markers <- head(c.markers[order(c.markers$avg_log2FC, decreasing = TRUE),],10)
  df <- bind_rows(df,c.markers)
  tiff(paste0("Cd8a/All_Cd8_cluster",as.character(x[i]),"_HM.tiff"),units="in", width=8, height=4.8, res=180)
  DefaultAssay(pbmc.integrated) <- "integrated"
  print(DoHeatmap(pbmc.integrated,features = rownames(df), draw.lines = TRUE) &
          scale_fill_gradientn(colours = oompaBase::redgreen(50)))
  #2
  tiff(paste0("Cd8a/All_Cd8_clusters_HM2.tiff",as.character(x[i]),"_HM.tiff"),units="in", width=8, height=8.8, res=180)
  DoHeatmap(subset(pbmc.integrated, cells =  WhichCells(pbmc.integrated, idents = x)),cells = WhichCells(pbmc.integrated, idents = x),features = df, draw.lines = TRUE,angle = 0) &
    scale_fill_gradientn(colours = oompaBase::redgreen(50)) & theme(text = element_text(size =8),legend.key.size = unit(.56, 'cm'))
  dev.off()
  c.markers <- cbind(c.markers,gene=rownames(c.markers))
  c_genes_names<-merge(c.markers, AnnotationDbi::select(org.Mm.eg.db, keys=as.character(c.markers$gene), columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL"), by.x="gene",by.y="SYMBOL") 
  
  
  
  # CP_res<-enrichGO(gene = as.vector(c_genes_names$ENTREZID),
  #                   OrgDb         = "org.Mm.eg.db",
  #                   ont           = "BP",
  #                   pAdjustMethod = "BH",
  #                   pvalueCutoff  = 0.01,
  #                 qvalueCutoff  = 0.2,
  #                   readable      = TRUE)
  
  
  
  tiff(paste0("Cd8a/All_Cd8_cluster",as.character(x[i]),"_GO.tiff"),units="in", width=7, height=8, res=180)
  print(clusterProfiler::dotplot(CP_res, showCategory = 20))
  dev.off()
  
  
  c.markers2 <- c.markers[order(c.markers$avg_log2FC, decreasing = TRUE),]
  #c.markers2 <- c.markers[names(sort(rowMeans(c.markers[,c(2,7,12,17)]), decreasing = TRUE)),]
  write.xlsx2(c.markers2, file= paste0("Cd8a/All_Cd8_cluster_",as.character(x[i]),".xlsx"))
}







