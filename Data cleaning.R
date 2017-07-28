###########################################
#		STAT 5010 Project EDA		#
# Predicting TE Preferences in V. fisheri #
###########################################

##Data cleaning
TEdata<-read.csv("TEmerged.csv", h=T)
head(TEdata)

TEmerged12<-read.csv("TEmerged12col.csv", h=T)
Headers<-which(TEmerged12[,1]=="Query id")
TEtest<-TEmerged12[-Headers,]
#Begin<-which(TEtest[,1]=="failed queries manual BLAST")
#Failed<-c(Begin-2,Begin-1, Begin)
#TEsave<-TEtest[-Failed,]
write.csv(TEtest,file="TEmerged12colblank.csv", row.names=F)

TE13<-read.csv("TEmerged13col.csv", h=T)
Headers<-which(TE13[,1]=="Query id")
TE13merged<-TE13[-Headers,]
write.csv(TE13merged, file="TEmerged13_nolabels.csv", row.names=F)

TE10<-read.csv("TEmerged10col.csv", h=T)
Headers<-which(TE10[,1]=="Query id")
Headers
TE10merged<-TE10[-Headers,]
write.csv(TE10merged, file="TEmerged10_nolabels.csv", row.names=F)

TEall<-read.csv(file.choose(),h=T)
Headers<-which(TEall[,1]=="Query.id")
TEall<-TEall[-Headers,]
write.csv(TEall, file="TEallmerged.csv", row.names=F)

TEallFinal<-read.csv(file.choose(),h=T)
head(TEallFinal)

#Chr 1 -- NC_006840.2
#Chr 2 -- NC_006841.2
#Plasm -- NC_006842.1

CN<-colnames(TEallFinal)

Ch1<-which(TEallFinal[,2]=="NC_006840.2")
Chrom1<-TEallFinal[Ch1,]
colnames(Chrom1)<-CN
write.csv(Chrom1, file="Ch1.csv", row.names=F)

Ch2<-which(TEallFinal[,2]=="NC_006841.2")
length(Ch2)
Chrom2<-TEallFinal[Ch2,]
colnames(Chrom2)<-CN
write.csv(Chrom2, file="Ch2.csv", row.names=F)

Plas<-which(TEallFinal[,2]=="NC_006842.1")
length(Plas)
Plasmid<-TEallFinal[Plas,]
colnames(Plasmid)<-CN
write.csv(Plasmid, file="Plamid.csv", row.names=F)




