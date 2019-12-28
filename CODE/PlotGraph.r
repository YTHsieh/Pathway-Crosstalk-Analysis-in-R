library(igraph)
datdir="D:/目標資料夾/"
links=read.delim(paste(datdir,"A組_crosstalk(full-unique).txt",sep=""),stringsAsFactors=F,header = T ,na.strings="NA",fill =T)
nodes=read.delim(paste(datdir,"A組_NODE.txt",sep=""),stringsAsFactors=F,header = T ,na.strings="NA",fill =T)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)    #把links和nodes整合入net檔(A組)

links=read.delim(paste(datdir,"B組_crosstalk(full-unique).txt",sep=""),stringsAsFactors=F,header = T ,na.strings="NA",fill =T)
nodes=read.delim(paste(datdir,"B組_NODE.txt",sep=""),stringsAsFactors=F,header = T ,na.strings="NA",fill =T)
net2 <- graph_from_data_frame(d=links, vertices=nodes, directed=F)    #把links和nodes整合入net檔(B組)


##-----參數調整測試--------------------------
plot(net)   ##嘗試plot出圖
E(net)   #net的edges


V(net)$size <- V(net)$exact.genes.NUM*0.02    #設定net的node大小及edge寬度
E(net)$width <- E(net)$JC*20

V(net2)$size <- V(net2)$exact.genes.NUM*0.02    #設定net2的node大小及edge寬度
E(net2)$width <- E(net2)$JC*20

###
l <- layout_in_circle(net)   ##設定固定layout的座標參數l
###

plot(net, edge.arrow.size=.2, edge.color="orange",
vertex.color="orange", vertex.frame.color="#ffffff",
vertex.label=V(net)$pathway, vertex.label.color="black", layout=l)     ##透過這組設定固定layout的座標，即可固定每次出圖的nodes位置(出A組)


plot(net2, edge.arrow.size=.2, edge.color="orange",
vertex.color="orange", vertex.frame.color="#ffffff",
vertex.label=V(net)$pathway, vertex.label.color="black", layout=l)     ##透過這組設定固定layout的座標，即可固定每次出圖的nodes位置(出B組)
##---------------------------------------------


##------並排出圖，左：A組；右：B組----------
par(mfrow=c(1,2), mar=c(0,0,0,0)) # plot two figures - 1 rows, 2 columns
plot(net, edge.arrow.size=.2, edge.color="orange",
vertex.color="orange", vertex.frame.color="#ffffff",
vertex.label=V(net)$pathway, vertex.label.color="black", layout=l) 
plot(net2, edge.arrow.size=.2, edge.color="orange",
vertex.color="orange", vertex.frame.color="#ffffff",
vertex.label=V(net)$pathway, vertex.label.color="black", layout=l) 
##---------------------------------------------

dev.off()     #關掉圖形
