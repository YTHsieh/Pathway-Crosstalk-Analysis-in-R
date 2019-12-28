##讀入pathway 資料
datdir <- "D:/目標資料夾目錄/"

target=read.delim(paste(datdir,"目標pathway.txt",sep = ""), sep = "\t",header = T) #讀入愈比較的目標pathway
pathway=read.delim(paste(datdir,"B組pathway內容.txt",sep=""),stringsAsFactors=F,header = T ,na.strings="NA",fill =T)	#讀入"B組"的「pathway內容」檔案
papaya=read.delim(paste(datdir,"A組pathway內容.txt",sep=""),stringsAsFactors=F,header = T ,na.strings="NA",fill =T)	#讀入"A組"的「pathway內容」檔案


##----此段為整理圖論的NODE資料------
for (i in 1:10){      
pathway[i,which(pathway[i,]=="")]="9"
}    ##空格補入9(A組和B組都要做一次)
pathway[10,2583]="9"   ##此格需額外補入，原因不明

NON=matrix(nrow=10,ncol=7)
NON[,1:6]=as.matrix(target[,1:6])
NON[,7]=as.matrix(pathway[,2584])     #重新計算pathway內的gene數(此部分運用了excel)
colnames(NON)=c("pathway.id","pathway","n.genes.in.pw","n.sex.related.genes.in.pw","p-value","p-value_BH_correction","exact.genes.NUM")
outfile=paste(datdir,"A組_NODE.txt", sep="");outfile	
write.table(NON,file=outfile, sep="\t",quote = F,row.names = F,col.names = T)	

##------------------------------------



##----此段為CROSSTALK分析-------------
A=list()     #創建list，讀入path(一個list一個path)
for (i in 1:10){
pathwayisNA=which((pathway[i,4:2199])=="9")                   ##"4:2199", 此數字為有gene的欄位數
A[[i]]=as.character(pathway[i,4:(pathwayisNA[1])])
}



BIG=matrix(nrow=100,ncol=6)     ##創建空matrix，以供後續分析

for (i in 1:10){

	for (k in 1:10){
	cat('now is handling the',i,k,'-th element (still',(10-i),'rows to go!!', '\n')
	INT=intersect(A[[i]],A[[k]])
	UNI=union(A[[i]],A[[k]])
	Ai.len=length(A[[i]])
	Ak.len=length(A[[k]])
	BIG[k+(i-1)*10,1:6]=c(pathway[i,1],pathway[k,1],length(INT),length(UNI),Ai.len,Ak.len)
	}
}
colnames(BIG)=c("Pathway1","Pathway2","Intersect","Union","Length(Ai)","Length(Ak)")    
###此步驟結束後，crosstalk分析做完(但未計算JC、OC值)

outfile=paste(datdir,"A組_crosstalk.txt", sep="");outfile	
write.table(BIG,file=outfile, sep="\t",quote = F,row.names = F,col.names = T)	#儲存階段檔。

#計算JC值
BIGG=matrix(nrow=100,ncol=8)
BIGG[,1:6]=BIG
colnames(BIGG)=c("Pathway1","Pathway2","Intersect","Union","Length(Ai)","Length(Ak)","JC","OC")
BIGG[,7]=as.numeric(BIG[,3])/as.numeric(BIG[,4])
BIGG[,7]=round(as.numeric(BIGG[,7]),digits = 3)   #取小數點下第三位

#計算OC值
for ( i in 1:100){		
		C=c(as.numeric(BIGG[i,5]),as.numeric(BIGG[i,6]))
		MIN=min(C)
		BIGG[i,8]=as.numeric(BIGG[i,3])/MIN	
		}
BIGG[,8]=round(as.numeric(BIGG[,8]),digits = 3)		
outfile=paste(datdir,"A組_crosstalk(full).txt", sep="");outfile	
write.table(BIGG,file=outfile, sep="\t",quote = F,row.names = F,col.names = T)	 #儲存階段檔。
			
			
#移除重複的crosstalk組合
BIGG=read.delim(paste(datdir,"papaya_crosstalk(full).txt",sep=""),stringsAsFactors=F,header = T ,na.strings="NA",fill =T)	
BIGGG=BIGG[-which(BIGG[,1]==BIGG[,2]),]
PAS12=paste(BIGGG[,1],BIGGG[,2],sep="")
PAS21=paste(BIGGG[,2],BIGGG[,1],sep="")
MINUS=BIGGG
MIno=c()
for (i in 1:90){
	cat('now is handling the',i,'-th element (still',(90-i),'rows to go!!', '\n')
	minusno=which(PAS12==PAS21[i])
	MIno=c(MIno,max(minusno,i))
	}
MIno2=unique(MIno);MIno2
BIGGGG=MINUS[-MIno2,]
outfile=paste(datdir,"A組_crosstalk(full-unique).txt", sep="");outfile	
write.table(BIGGGG,file=outfile, sep="\t",quote = F,row.names = F,col.names = T)	
			
#完成全部參數計算。

##---------------------------------
