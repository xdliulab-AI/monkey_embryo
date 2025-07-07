#3D gradient construction
fun_color_range <- colorRampPalette(c("blue","#e0e0e0", "#b2182b"))
cs9metadata <- read.csv("CS10_SPATEO_align_metadata.csv",header = TRUE,row.names = 1)
unique((cs9metadata$celltype))
cs9metadata_embryo <- cs9metadata[cs9metadata$celltype%in%c("Brain and spinal cord","Heart","Gut Tube","Paraxial Mesoderm","Endothelial cells","Lateral Plate Mesoderm","Allantois","Caudal Mesoderm","Notochord"),]

plot3d(x=cs9metadata_embryo$x,y=cs9metadata_embryo$y,z=cs9metadata_embryo$z)   

euclidean_distance <- function(point1, point2) {
  return(sqrt(sum((point2 - point1)^2)))
}
ini_p <- c(50,100,6)
#allpoint <- cs9metadata_embryo[1:1000,c("x","y","z")]
allpoint <- cs9metadata_embryo[,c("x","y","z")]
rownames(allpoint) <- paste0(allpoint$x,"_",allpoint$y,"_",allpoint$z)

pointstmp <- allpoint[,c("x","y","z")]
points <- as.matrix(allpoint[,c("x","y","z")])
#percent <- nrow(pointstmp) / total_rows
#print(paste0("percent-",percent*100,"%"))
print(length(pointstmp$x))
reference_point <- matrix(ini_p, nrow = 1, ncol = 3, byrow = TRUE)
distances <- apply(points, 1, function(point) {
  euclidean_distance(reference_point, point)
})
pointstmp$distances <- distances
pointstmp <- pointstmp[order(pointstmp$distances),]
dim(pointstmp)

pointstmp$color <- fun_color_range(33887) 
library(rgl)
plot3d(x=pointstmp$x,y=pointstmp$y,z=pointstmp$z,col = pointstmp$color,zlim = c(0,30))    

rankn <- 200
# 创建示例数据框
df <- pointstmp
df <- df[order(df$distances),]

n <- nrow(df)
size_per_portion <- ceiling(n / rankn)

portion_number <- 1
df$Aportion_id <- NA

for (i in 1:n) {
  df$Aportion_id[i] <- portion_number
  if (i %% size_per_portion == 0 && portion_number < rankn) {
    portion_number <- portion_number + 1
  }
}

library(plyr)
colors <- fun_color_range(200) 
length(unique(df$Aportion_id))
df$color=as.vector(mapvalues(df$Aportion_id,from = unique(df$Aportion_id),to=colors))

head(df)
plot3d(x=df$x,y=df$y,z=df$z,col = df$color,zlim = c(0,30))

df$Aportion_names <- paste0("AP_",df$Aportion_id)

#计算每一个AP的中位数点作为中心线
xmean <- aggregate(x ~ Aportion_names, data = df, FUN = median)
ymean <- aggregate(y ~ Aportion_names, data = df, FUN = median)
zmean <- aggregate(z ~ Aportion_names, data = df, FUN = median)
center_curve1 <- merge(xmean,ymean,by = "Aportion_names")
center_curve <- merge(center_curve1,zmean,by = "Aportion_names")
# 对中心曲线进行平滑处理（这里使用 loess 回归进行简单平滑）
#smooth_curve <- predict(loess(y ~ x, data = center_curve))
#smooth_center_curve <- data.frame(x = xmean$x, y = smooth_curve)
library(ggplot2)

plotlist <- list()
for (i in unique((df$z))) {
  plotlist[[as.character(i)]] <- ggplot()+geom_point(data =df[df$z==i,], aes(x,y))+geom_point(data = center_curve,aes(x,y),color = "red")+coord_fixed()+theme_classic()
}
patchwork::wrap_plots(plots = plotlist)
p1 <- ggplot()+geom_point(data =df[df$z==0,], aes(x,y))+geom_point(data = center_curve,aes(x,y),color = "red")+coord_fixed()+theme_classic()
p1
rownames(center_curve) <- center_curve$Aportion_names
center_curve$z <- 5.4
# 计算每个点到中心曲线的距离
distances <- apply(df[,c("Aportion_names","x","y")], 1, function(row) {
  Apid <-  as.character(row[1])
  point <- as.numeric(row[2:3])
  nearest_point_on_curve <- as.numeric(center_curve[Apid,c("x","y")])
  #print(Apid)
  #print(point)
  mark <- (nearest_point_on_curve[1] - point[1])/abs(nearest_point_on_curve[1] - point[1])
  dd <- euclidean_distance(point, nearest_point_on_curve)
  return(dd*mark)
})
distances <- as.data.frame(distances)
distances$x <- df[rownames(distances),"x"]
distances$y <- df[rownames(distances),"y"]
distances$z <- df[rownames(distances),"z"]
distances[is.na(distances),"distances"] <- 0
distances$Aportion_names <- df[rownames(distances),"Aportion_names"]
data <- distances
library(dplyr)

# 按组进行排序
data_sorted <- data %>%
  group_by(Aportion_names) %>%
  arrange(distances)

# 计算每个组的总数量
group_sizes <- data_sorted %>%
  group_by(Aportion_names) %>%
  summarise(size = n())
#group_sizes$size <- 100

# 计算每个组中每份的大小
group_divisions <- group_sizes %>%
  mutate(step_size = size / 100)

#group_divisions <- 100
# 为每个数据点分配所属的份数
data_divided <- data_sorted %>%
  left_join(group_divisions, by = "Aportion_names") %>%
  mutate(division = ceiling(row_number() / step_size))

# 查看结果
data_divided <- as.data.frame(data_divided)
rownames(data_divided) <- paste0(data_divided$x,"_",data_divided$y,"_",data_divided$z)
df$Ddistance <- data_divided[rownames(df),"distances"]
df$Dportion_id <- data_divided[rownames(df),"division"]


fun_color_range <- colorRampPalette(c("blue","#e0e0e0", "#b2182b"))
length(sort(unique(df$Dportion_id)))
df$color=as.vector(mapvalues(df$Dportion_id,from = sort(unique(df$Dportion_id)),to=fun_color_range(101)))

plot3d(x=df$x,y=df$y,z=df$z,col = df$color,zlim = c(0,10))

cs9metadata$index <- paste0(cs9metadata$x,"_",cs9metadata$y,"_",cs9metadata$z)
cs9metadata$cellLable <- rownames(cs9metadata)
df$index <- rownames(df)

df_result <- merge(df,cs9metadata,by = "index")
