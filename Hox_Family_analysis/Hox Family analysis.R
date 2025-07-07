#Hox Family analysis
library(ggplot2)
library(Matrix)
library(Seurat)
library(data.table)

mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 12),
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 11),
                 axis.text.x = element_text(color = "black"),
                 axis.text.y = element_text(color = "black")) #自定义主题

HOXC1 <- list("C1"=c("HOXA1","HOXB1","HOXD1"))
HOXC2 <- list("C2"=c("HOXA2","HOXA3","HOXB2","HOXB3","HOXD3"))
HOXC3 <- list("C3"=c("HOXA4","HOXB4","HOXC4","HOXD4"))
HOXC4 <- list("C4"=c("HOXA5","HOXB5","HOXC5"))
HOXC5 <- list("C5"=c("HOXA6","HOXB6","HOXC6"))
HOXC6 <- list("C6"=c("HOXA7","HOXB7"))
HOXC7 <- list("C7"=c("HOXB8","HOXC8","HOXD8"))
HOXC8 <- list("C8"=c("HOXA9","HOXB9","HOXC9","HOXD9","HOXA10","HOXC10","HOXD10","HOXA11","HOXC11","HOXD11","HOXC12","HOXD12","HOXA13","HOXB13","HOXC13","HOXD13"))

obj <- readRDS("Ddistance.rds")

obj <- AddModuleScore(obj,
                      features = c(HOXC1,HOXC2,HOXC3,HOXC4,HOXC5,HOXC6,HOXC7,HOXC8),
                      ctrl = 100)

ng <- HOXC1[["C1"]]
length(ng)
metadatatmp <- cbind(obj@meta.data[,c("x","y","Aportion_id","Aportion_names")],as.data.frame(t(as.data.frame(obj@assays$RNA@counts[ng,]))))
names(metadatatmp) <- c("x","y","Aportion_id","Aportion_names",ng)
result <- list()
for (gene in ng){
  print(gene)
  Sys.sleep(0.1)
  df <- metadatatmp
  result[[gene]] <- df %>%
    group_by(Aportion_names) %>%
    summarise(gt_zero = sum(get(gene) > 0),
              eq_zero = sum(get(gene) == 0))
  result[[gene]]$portion <- result[[gene]]$gt_zero/(result[[gene]]$eq_zero+result[[gene]]$gt_zero)
  result[[gene]]$Adistancerank <- as.integer(sub(".*_","",result[[gene]]$Aportion_names))
  names(result[[gene]]) <- c("Aportion_names","gt_zero","eq_zero",'portion',"Adistancerank")
  result[[gene]] <- as.data.frame(result[[gene]])
  rownames(result[[gene]]) <- result[[gene]]$Aportion_names
  result[[gene]] <- result[[gene]][,c("gt_zero","eq_zero",'portion')]
  result[[gene]]$gene <- gene
  #result[[gene]]$gene <- gene
  print(which(ng == gene))
}
allresult <- do.call(rbind,result)
allresult$Adistancerank <- as.integer(sub(".*_","",rownames(allresult)))

p1 <- ggplot(allresult,aes(Adistancerank,portion,color = gene))+geom_point(size=0.5)+geom_smooth(aes(fill = gene),method = 'loess',se=TRUE,
                                                                                                 color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))
