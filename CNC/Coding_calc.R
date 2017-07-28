###############################################
# Coding and Non-coding regions of the genome #
###############################################
Ch1Gene<-read.table("Ch1Gene.txt",sep="\t", h=T, quote="")
Ch2Gene<-read.table("Ch2Gene.txt",sep="\t", h=T, quote="")
PlsmGene<-read.table("PlsmGene.txt",sep="\t", h=T, quote="")

#############################################
# Function that takes in a NCBI data frame  #
# and calculates the length of genes and    #
# the distance between each gene; returns   #
# a new data frame with the results cbinded #
#############################################
LengthDiff<-function(Data){
	Data<-Data[order(Data[,13]),]
	row.names(Data)<-1:nrow(Data)

	gene_length<-Data[,14]-Data[,13] #Length of the gene
	
	gene_dist<-vector(length=nrow(Data))
	i=1
	while(i<nrow(Data)){ #Calculates the distance between genes; can be negative if overlap
		gene_dist[i]<-Data[i+1,13]-Data[i,14]
		i=i+1
	} #End of distance while-loop
Result<-cbind(Data, gene_length, gene_dist)
return(Result)
}#End of function


#Calculating for each genetic element
C1<-LengthDiff(Ch1Gene)
C2<-LengthDiff(Ch2Gene)
P<-LengthDiff(PlsmGene)

#################################################
# Function to take in the lengths and distances #
# of genes to calculate a genome vector         #
#################################################
CNC<-function(Data){
	Genome<-vector()
	
	for(i in 1:nrow(Data)){
		if(Data[i,20]<0){ #If a gene overlaps, shorten the first gene
			coding<-rep(1, times=(Data[i,19]-Data[i,20]))
			non<-vector(length=0)
		}#End if

		else{ #If the genes do not overlap
			coding<-rep(1, times=Data[i,19])
			non<-rep(0, times=Data[i,20])
		}#End else
	Genome<-c(Genome, coding, non)
	}#End for-loop

return(Genome)
}#End function


#Calculating coding vs non-coding genome vectors and proportions
C1Genome<-CNC(C1)
(sum(C1Genome))/(length(C1Genome)) #0.8900263
write.table(C1Genome, file="Ch1_CodingSeq.txt",row.names=F)

C2Genome<-CNC(C2)
(sum(C2Genome))/(length(C2Genome)) #0.8816424
write.table(C2Genome, file="Ch2_CodingSeq.txt", row.names=F)

PGenome<-CNC(P)
(sum(PGenome))/(length(PGenome)) #0.87646
write.table(PGenome, file="Plasmid_CodingSeq.txt", row.names=F)

