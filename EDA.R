#################################
# EDA for STAT Capstone Project #
#################################
library(ggplot2)
library(psych)

##################
#  Chromosome 1  #
##################

Ch1<-read.csv("Ch1.csv", h=T)
Ch1<-read.csv(file.choose(),h=T)
sort(table(Ch1[,9])) #Most freq - VF_0328 with 154 hits
Length1<-data.frame(Ch1[,4]-Ch1[,3])

#Alginment starting positions
png("AlignmentStarting_Ch1.png", width=1600, height=900)
	ggplot(Ch1, aes(Ch1[,5])) + geom_histogram() +
	ggtitle("Distribution of alignment positions in chromosome 1") + 
	labs(x="Alignment Positions", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Ch1[,5]))

#Frequency of the genes the TE is inserted in
png("GeneFreq_Ch1.png", height=900, width=1600)
ggplot(Ch1, aes(Ch1[,9])) + geom_bar() +
	ggtitle("Frequency of insertion genes in chromsome 1") +
	labs(x="Gene", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Ch1[,9]))

#Distribution of alignment lengths
png("AlignmentDist_Ch1.png", height=900, width=1600)
ggplot(Length1, aes(Length1[,1])) + geom_histogram() +
	ggtitle("Distribution of alignment lengths in chromsome 1") +
	labs(x="Alignment Length", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
describe(Length1)
sort(table(Length1))

#Frequency of the strand
png("StrandFreq_Ch1.png", height=900, width=1600)
	barplot(table(Ch1[,8]), main="Frequency of strands")
dev.off()
table(Ch1[,8])


##################################
# Figures for report - no titles #
##################################
#Alginment starting positions
png("Report_TEinsertions_Ch1.png")
	ggplot(Ch1, aes(Ch1[,5])) + geom_histogram(binwidth=75000) +
	labs(x="Position on Chromosome (bp)", y="Number of Insertions") + theme(plot.title=element_text(size=32, face="bold"),
	axis.text=element_text(size=16), axis.title=element_text(size=26, face="bold", hjust=.5))
dev.off()
sort(table(Ch1[,5]))

#Frequency of the genes the TE is inserted in
png("Report_GeneFreq_Ch1.png", height=900, width=1600)
ggplot(Ch1, aes(Ch1[,9])) + geom_bar() +
#	ggtitle("Frequency of insertion genes in chromsome 1") +
	labs(x="Gene", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Ch1[,9]))

#Distribution of alignment lengths
png("Report_AlignmentDist_Ch1.png", height=900, width=1600)
ggplot(Length1, aes(Length1[,1])) + geom_histogram() +
#	ggtitle("Distribution of alignment lengths in chromsome 1") +
	labs(x="Alignment Length", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
describe(Length1)
sort(table(Length1))

#Frequency of the strand
png("Report_StrandFreq_Ch1.png", height=900, width=1600)
	barplot(table(Ch1[,8])) #, main="Frequency of strands")
dev.off()
table(Ch1[,8])




################################################################################
# This function takes in a data frame assumming exact column order of original #
# and calculates the frequency of the insertion of the TE in the gene for each #
# gene that has at least 2 insertions in the entire chromosome/plasmid         #
# A new folder containing the output files is made					 #
# 	Note: be sure to use "" when assigning the name of the directory 		 #
#		using the variable Chromdrir                                       #
################################################################################
InsertionGene<-function(Data, Chromdir){ 
library(ggplot2)
library(gridExtra)
dir.create(Chromdir)
Levels<-which(table(Data[,9])>=2) #Vector of all the genes that have >= 2 insertions

	for(i in 1:length(Levels)){ #For each gene
		Gene<-which(Data[,9]==names(Levels[i])) 
		DataGeneSub<-Data[Gene,] #Retrieve rows corresponding to the gene
		InsertPos<-data.frame(DataGeneSub[,5]-DataGeneSub[,11]) #Calculate insertion position
	
		#Plotting results
		G<-ggplot(InsertPos, aes(x=InsertPos[,1])) + geom_histogram() +
			ggtitle(paste0("Frequency of TE Insertions in ", names(Levels[i]))) + 
			labs(x=paste0("Position in ", names(Levels[i]), "gene"), y="Frequency") + 
			theme(plot.title=element_text(size=44, face="bold"),
			axis.text=element_text(size=40), 
			axis.title=element_text(size=40,face="bold", hjust=.5))

		png(file=paste0(Chromdir, "/", names(Levels[i]), "_FreqPlot.png"), width=1600, height=900)
		grid.arrange(G,nrow=1,ncol=1)
		dev.off()		


		ReportPlot<-ggplot(InsertPos, aes(x=InsertPos[,1])) + geom_histogram() +
			#ggtitle(paste0("Frequency of TE Insertions in ", names(Levels[i]))) + 
			labs(x=paste0("Position in ", names(Levels[i]), "gene"), y="Frequency") + 
			theme(plot.title=element_text(size=44, face="bold"),
			axis.text=element_text(size=40), 
			axis.title=element_text(size=40,face="bold", hjust=.5))

		png(file=paste0(Chromdir, "/Report_", names(Levels[i]), "_FreqPlot.png"), width=1600, height=900)
		grid.arrange(ReportPlot,nrow=1,ncol=1)
		dev.off()
	
		#Saving results
		GeneTableFreq<-(sort(table(InsertPos)))
		write.csv(GeneTableFreq, file=paste0(Chromdir, "/", names(Levels[i]),"_Table.csv"), row.names=F)				
	}#End of levels for loop
	
	#Outputs the file that has the exact insertion location in the gene and the frequencies of insertions
write.csv(t(Levels), file=paste0(Chromdir, "/TotalGenes_Insertions.csv"), row.names=F)
} #End of the function

InsertionGene(Plasmid, "Plasmid_InsertionFreq")
InsertionGene(Ch1, "Ch1_InsertionFreq")
InsertionGene(Ch2, "Ch2_InsertionFreq")


#####################
# Percent function #
#####################
GenePercent<-function(Data, directory){
library(reshape2)
Percents<-vector()

Levels<-which(table(Data[,9])>=1) #Vector of all the genes that have >= 1 insertions

	for(i in 1:length(Levels)){ #For each gene
		Gene<-which(Data[,9]==names(Levels[i])) 
		DataGeneSub<-Data[Gene,] #Retrieve rows corresponding to the gene
		InsertPos<-data.frame(DataGeneSub[,5]-DataGeneSub[,11])
		GeneLength<-DataGeneSub[,12]-DataGeneSub[,11]
	
		InsertPercent<-InsertPos/GeneLength
	
		Melted<-melt(InsertPercent)
		Percents<-c(Percents, Melted$value) 
	}
png("Insert_GenePercent_Ch1.png")
hist(Percents, breaks="Scott", xlab="Gene Length Percentage", main="", cex.axis=1.4)
dev.off()

}

S<-GenePercent(Ch1)
S
##################
#  Chromosome 2  #
##################

Ch2<-read.csv("Ch2.csv", h=T)
Ch2<-read.csv(file.choose(),h=T)
sort(table(Ch2[,9])) #Most freq - VF_A0451 with 74 hits
Length2<-data.frame(Ch2[,4]-Ch2[,3])
describe(Length2)
sort(table(Length2))

#TE insertion dist across ch2
Ch2Length<-1:1330333
png("InsertionDist_Ch2.png")
	ggplot(Ch2, aes(Ch2[,5])) + geom_histogram()
dev.off()
hist(Ch2[,5])

#Alginment starting positions
png("AlignmentStarting_Ch2.png", width=1600, height=900)
	ggplot(Ch2, aes(Ch2[,5])) + geom_histogram() +
	ggtitle("Distribution of alignment positions in chromosome 2") + 
	labs(x="Alignment Positions", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Ch2[,5]))

#Frequency of the genes the TE is inserted in
png("GeneFreq_Ch2.png", height=900, width=1600)
ggplot(Ch2, aes(Ch2[,9])) + geom_bar() +
	ggtitle("Frequency of insertion genes in chromsome 2") +
	labs(x="Gene", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Ch2[,9]))

#Distribution of alignment lengths
png("AlignmentDist_Ch2.png", height=900, width=1600)
ggplot(Length2, aes(Length2[,1])) + geom_histogram() +
	ggtitle("Distribution of alignment lengths in chromsome 2") +
	labs(x="Alignment Length", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()

#Frequency of the strand
png("StrandFreq_Ch2.png", height=900, width=1600)
	barplot(table(Ch2[,8]), main="Frequency of strands")
dev.off()
table(Ch2[,8])


#################################
# Figures for report - no title #
#################################
#Alginment starting positions
png("Report_TEinsertions_Ch2.png")
	ggplot(Ch2, aes(Ch2[,5])) + geom_histogram(binwidth=50000) +
	labs(x="Position on Chromosome (bp)", y="Number of Insertions") + theme(plot.title=element_text(size=32, face="bold"),
	axis.text=element_text(size=16), axis.title=element_text(size=26, face="bold", hjust=.5))
dev.off()
sort(table(Ch2[,5]))

#Frequency of the genes the TE is inserted in
png("Report_GeneFreq_Ch2.png", height=900, width=1600)
ggplot(Ch2, aes(Ch2[,9])) + geom_bar() +
#	ggtitle("Frequency of insertion genes in chromsome 2") +
	labs(x="Gene", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Ch2[,9]))

#Distribution of alignment lengths
png("Report_AlignmentDist_Ch2.png", height=900, width=1600)
ggplot(Length2, aes(Length2[,1])) + geom_histogram() +
#	ggtitle("Distribution of alignment lengths in chromsome 2") +
	labs(x="Alignment Length", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()

#Frequency of the strand
png("Report_StrandFreq_Ch2.png", height=900, width=1600)
	barplot(table(Ch2[,8])) #, main="Frequency of strands")
dev.off()
table(Ch2[,8])


#############
#  Plasmid  #
#############

Plasmid<-read.csv("Plasmid.csv", h=T)
Plasmid<-read.csv(file.choose(), h=T)
sort(table(Plasmid[,9])) #Most freq - VF_B0050 with 46 hits
LengthP<-data.frame(Plasmid[,4]-Plasmid[,3])
describe(LengthP)
sort(table(LengthP))

#Alginment starting positions
png("AlignmentStarting_Plasmid.png", width=1600, height=900)
	ggplot(Plasmid, aes(Plasmid[,5])) + geom_histogram() +
	ggtitle("Distribution of alignment positions in the plasmid") + 
	labs(x="Alignment Positions", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Plasmid[,5]))

#Frequency of the genes the TE is inserted in
png("GeneFreq_Plasmid.png", height=900, width=1600)
ggplot(Plasmid, aes(Plasmid[,9])) + geom_bar() +
	ggtitle("Frequency of insertion genes in the plasmid") +
	labs(x="Gene", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Plasmid[,9]))

#Distribution of alignment lengths
png("AlignmentDist_Plasmid.png", height=900, width=1600)
ggplot(LengthP, aes(LengthP[,1])) + geom_histogram() +
	ggtitle("Distribution of alignment lengths in the plasmid") +
	labs(x="Alignment Length", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(LengthP[,1]))

#Frequency of the strand
png("StrandFreq_Plasmid.png", height=900, width=1600)
	barplot(table(Plasmid[,8]), main="Frequency of strands")
dev.off()
table(Plasmid[,8])


#################################
# Figures for report - no title #
#################################
#Alginment starting positions
png("Report_TEinsertions_Plasmid.png")
	ggplot(Plasmid, aes(Plasmid[,5])) + geom_histogram(binwidth=2000) +
#	ggtitle("Distribution of alignment positions in the plasmid") + 
	labs(x="Position on Plasmid (bp)", y="Number of Insertions") + theme(plot.title=element_text(size=32, face="bold"),
	axis.text=element_text(size=16), axis.title=element_text(size=26, face="bold", hjust=.5))
dev.off()
sort(table(Plasmid[,5]))

#Frequency of the genes the TE is inserted in
png("Report_GeneFreq_Plasmid.png", height=900, width=1600)
ggplot(Plasmid, aes(Plasmid[,9])) + geom_bar() +
#	ggtitle("Frequency of insertion genes in the plasmid") +
	labs(x="Gene", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(Plasmid[,9]))

#Distribution of alignment lengths
png("Report_AlignmentDist_Plasmid.png", height=900, width=1600)
ggplot(LengthP, aes(LengthP[,1])) + geom_histogram() +
#	ggtitle("Distribution of alignment lengths in the plasmid") +
	labs(x="Alignment Length", y="Frequency") + theme(plot.title=element_text(size=44, face="bold"),
	axis.text=element_text(size=40), axis.title=element_text(size=40, face="bold", hjust=.5))
dev.off()
sort(table(LengthP[,1]))

#Frequency of the strand
png("Report_StrandFreq_Plasmid.png", height=900, width=1600)
	barplot(table(Plasmid[,8])) #, main="Frequency of strands")
dev.off()
table(Plasmid[,8])




