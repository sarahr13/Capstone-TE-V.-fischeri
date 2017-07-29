PGC<-read.table("Plasmid_GCcontentSeq.txt", h=TRUE)
PCoding<-read.table("Plasmid_CodingNoncodingSeq.txt", h=TRUE)
Plasmid<-read.csv("Plasmid.csv",h=TRUE)

Ch2GC<-read.table("Ch2_GCcontentSeq.txt", h=TRUE)
Ch2Coding<-read.table("Ch2_CodingNoncodingSeq.txt", h=TRUE)
Ch2<-read.csv("Ch2.csv", h=TRUE)

Ch1GC<-read.table("Ch1_GCcontentSeq.txt", h=TRUE)
Ch1Coding<-read.table("Ch1_CodingNoncodingSeq.txt", h=TRUE)
Ch1<-read.csv("Ch1.csv",h=TRUE)

###################################################################################################
# Function calculates binwidth data - # insertions, GC-content prop, & coding vs. non-coding prop #
# Input: data frame containing the insertion data, and vectors for gc-content, and coding         #
# Output: list of data frames containing # insertions, prop of gc & coding content in each bin    #
# Note: The last bin is thrown out if it is not the exact binwidth                                #
###################################################################################################
bin_counter<- function(data, binwidths, gcContent, CodingContent){
	BinResults<-list() #Stores final results for all binwidths
	for(k in 1:length(binwidths)){ #For each binwidth
		binwidth<-binwidths[k]
	
		#Determines number of bins needed and assigns which bin each insertion goes into
		Bins<-vector(length=(ceiling(max(data[,6])/binwidth)))
		InsertBins<-ceiling(data[,5]/binwidth)
		count0=table(factor(InsertBins, levels=c(1:length(Bins))))

	#Function to split up the vector into correct binwidths
		splitAt2 <- function(x, pos) { 
		    out <- list()
		    pos2 <- c(1, pos, length(x)+1)
		    for (i in seq_along(pos2[-1])){
	       	 out[[i]] <- x[pos2[i]:(pos2[i+1]-1)]
          	    }#End of for loop
	  	return(out)
	}#End of split function

	#Calculate bin break locations
		BinBreaks<-vector()
		i=1
		while((i+binwidth)<length(gcContent[,1])){
			Pos<-i+binwidth
			BinBreaks<-c(BinBreaks, Pos)
			i=i+binwidth
		}#End of calculating bin breaks

	#GC-content
		GCsplit<-splitAt2(gcContent[,1], BinBreaks)
		GCprop<-unlist(lapply(GCsplit, mean))

	#Coding percentage
		CodingSplit<-splitAt2(CodingContent[,1], BinBreaks)
		CodingProp<-unlist(lapply(CodingSplit, mean))
	
	#Saving results from specified binwidth
	OneBin<-cbind(InsertionFreq=count0, GCprop, CodingProp)
	BinResults[[paste0("Binwidth",binwidth)]]<-OneBin	

	}#End of loop over different bin widths

return(BinResults)
}#End of entire bin calculator function

Ch1Bins<-c(1000, 5000, 10000)
TestCh1<-bin_counter(Ch1, Ch1Bins, Ch1GC, Ch1Coding)

Ch2Bins<-c(1000, 5000, 10000)
TestCh2<-bin_counter(Ch2, Ch2Bins, Ch2GC, Ch2Coding)

PlasmidBins<-c(750, 1000, 1500)
TestPlasmid<-bin_counter(Plasmid, PlasmidBins, PGC, PCoding)

###############################################
# Creating a simulated genome with insertions #
# Simulated insertion locations from a runif  #
# with the max and min as the beginning & end #
# of the chromosome				                    #
###############################################
Simbin_counter<- function(data, binwidths, Length){
	BinResults<-list() #Stores final results for all binwidths
	for(k in 1:length(binwidths)){ #For each binwidth
		binwidth<-binwidths[k]
		
		Bins<-vector(length=(ceiling(max(Length)/binwidth)))
		InsertBins<-ceiling(data/binwidth)
		InsertBins
		count0=table(factor(InsertBins, levels=c(1:length(Bins))))
		
		#Saving results from specified binwidth
		OneBin<-cbind(InsertionFreq=count0)
		BinResults[[paste0("Binwidth",binwidth)]]<-OneBin	

	}#End of for-k-loop
return(BinResults)
} #End of function

#############################################################################################
# Function that replicates the simulation a specified number of times 				              #
# element: comes from the bin_calculator function -- ONLY input the bin width you want!     #
# elem_length: length of the chromosome/plasmid in bp 						                          #
# binwidth: the specified bin width; should be equal to the binwidth of the actual provided #
# numreps: how many replicates										                                          #
# numActual: number of insertions in the chromosome/plasmid						                      #
# plotTitle: name of the title when it's output as a png						                        #
#############################################################################################
new_sim<-function(element, elem_length, binwidth, numreps=1000, NumActual, plotTitle, CSVname){
	SignifBins<-NULL
	Result<-list()
	for(j in 1:3){ #Go through each bin width
		dataf<-NULL
		Counts<-NULL
		Locations<-vector()
		Counts<-element[[j]][,1] #Changing bin widths
	
		for(i in 1:numreps){ #Create n binned replicates
		  seq<-1:elem_length
		  Sample<-sample(seq, size=NumActual, replace=TRUE)
		  simrep<-Simbin_counter(Sample, binwidth[j], elem_length)
		  Inserts<-simrep[[1]][,1]
		  dataf<-cbind(dataf, Inserts)
		}
	extremevalues<-NULL
	means<-NULL
	
	colMax <- function(data) sapply(data, max, na.rm = TRUE)
	extremevalues<-colMax(dataf)
	
	mu<-mean(dataf)
	cutoff<-mean(head(sort(extremevalues, decreasing=TRUE), floor(.05*length(Counts))))
	
	Bins<-length(which(Counts>cutoff)) #Calculating how many bins at that bin width are significant
	SignifBins[j]<-Bins
	Locations<-data.frame(binwidth[j]*(which(Counts>cutoff)))

	Labels<-seq(from=1, to=length(Counts), by=1)
	png(paste0(plotTitle, "_Binwidth",binwidth[j], ".png"))
	barplot(Counts, xlab="Bins",ylab="Insertions", cex.axis=1.4, names.arg=Labels)
	abline(h=cutoff,col='blue')

	abline(h=mu,col='red')	
	dev.off()

	#Saving the result for each bin width
		Result[paste0("Binwidth", binwidth[j], "_", SignifBins[j], "SigBins")]<-Locations
	}#End of loop over each bin

#Recording the locations of interest - add the bin width to the location to get the total bin width range
library(qpcR)
df<-qpcR:::cbind.na(Result[[1]],Result[[2]], Result[[3]])
colnames(df)<-c(paste0(names(Result)[1]),
			paste0(names(Result)[2]),
			paste0(names(Result)[3]))
write.csv(df, file=paste0(CSVname, "_SigBin_Locations.csv"), na="", row.names=FALSE)
return(df)
}##End of function

S<-new_sim(TestCh1, nrow(Ch1Coding), Ch1Bins, 5000, nrow(Ch1), "Simulated_Ch1", "Ch1")
S2<-new_sim(TestCh2, nrow(Ch2Coding), Ch2Bins, 5000, nrow(Ch2), "Simulated_Ch2", "Ch2")
S3<-new_sim(TestPlasmid , nrow(PCoding), PlasmidBins, 5000, 207, "Simulated_Plasmid", "Plasmid")


##############################################################################
png("CodingContentCh1.png")
boxplot(TestCh1[[1]][,3], TestCh1[[2]][,3], TestCh1[[3]][,3], names=c("1000", "5000", "10000"), xlab="Binwidth (bp)", ylab="Coding percentage", cex.lab=1.4)
dev.off()

png("CodingContentCh2.png")
boxplot(TestCh2[[1]][,3], TestCh2[[2]][,3], TestCh2[[3]][,3], names=c("1000", "5000", "10000"), xlab="Binwidth (bp)", ylab="Coding percentage", cex.lab=1.4)
dev.off()

#GC-content comparison between Ch1 and Ch2
par(mfrow=c(1,2))
boxplot(TestCh1[[1]][,2], TestCh1[[2]][,2], TestCh1[[3]][,2], names=c("1000", "5000", "10000"), xlab="Binwidth (bp)", ylab="Coding percentage", cex.lab=1.4)
boxplot(TestCh2[[1]][,2], TestCh2[[2]][,2], TestCh2[[3]][,2], names=c("1000", "5000", "10000"), xlab="Binwidth (bp)", ylab="Coding percentage", cex.lab=1.4)

#Coding percentage comparison between the chromosomes
par(mfrow=c(1,2))
boxplot(TestCh1[[1]][,3], TestCh1[[2]][,3], TestCh1[[3]][,3], names=c("1000", "5000", "10000"), xlab="Binwidth (bp)", ylab="Coding percentage", cex.lab=1.4)
boxplot(TestCh2[[1]][,3], TestCh2[[2]][,3], TestCh2[[3]][,3], names=c("1000", "5000", "10000"), xlab="Binwidth (bp)", ylab="Coding percentage", cex.lab=1.4)


###################################################################################
# Function that performs NEGATIVE BINOMIAL regression on the number of insertions #
# with the percentages of GC-content and coding content proportion in each bin as #
# the predictor variables. Returns the p-values                                   #
###################################################################################
NBreg<-function(DataList){
	library(MASS)
	Result<-list()
	for(i in 1:length(DataList)){
		InsertionY<-DataList[[i]][,1]
		gcX<-DataList[[i]][,2]
		codingX<-DataList[[i]][,3]

		#Regression w/ no interaction term
		#bothreg<-glm.nb(InsertionY~gcX+codingX)
		#bothRegSum<-summary(glm.nb(InsertionY~gcX+codingX))
		
		#Regression including interaction term
		bothInt<-glm.nb(InsertionY~gcX*codingX)
		bothIntSum<-summary(glm.nb(InsertionY~gcX*codingX))

		#IndResult<-list("noI"=bothreg, "interaction"=bothInt, "noInteractionSummary"=bothRegSum, "interactionSummary"=bothIntSum)
		IndResult<-list("interaction"=bothInt, "interactionSummary"=bothIntSum)
		Result[[i]]<-IndResult
	}#End of for each bin loop
return(Result)
}


T1<-NBreg(TestCh1)
T1

T2<-NBreg(TestCh2)
T2

T<-NBreg(TestPlasmid)
T


par(mfrow=c(1,2))
plot(log(predict(T[[1]][[2]], type="response")),T[[1]][[2]]$residuals)
plot(predict(T[[1]][[2]]), T[[1]][[2]]$residuals)

#Testing to see which model showed whether or not the predicted values 
#were significant in predicting the actual counts 
#Change the first pair of brackets to change the binwidth
summary(lm(TestPlasmid[[1]][,1]~predict(T[[1]][[1]])))
summary(lm(TestPlasmid[[1]][,1]~predict(T[[1]][[2]])))

##########################################
# 	Cross-validation on the NBR          #
# Input should come from the df that is  #
# output by the bins_calculator function #
##########################################
CV<-function(DataList, k=5, r=10){
library(MASS)
library(cvTools)
cvR<-list()
for(i in 1:length(DataList)){
	InsertionY<-DataList[[i]][,1]
	gcX<-DataList[[i]][,2]
	codingX<-DataList[[i]][,3]
	
	#Regression w/ both as predictors
	bothreg<-glm.nb(InsertionY~gcX+codingX, maxit=100)		
			
	#Regression w/ both as predictors, including interaction term
	bothInt<-glm.nb(InsertionY~gcX*codingX, maxit=100)

	noIntCV<-cvFit(bothreg, data=as.data.frame(DataList[[i]]), y=DataList[[i]][,1], cost=rtmspe,
			K=k, R=r, costArgs=list(trim=0.1))
	
	IntCV<-cvFit(bothInt, data=as.data.frame(DataList[[i]]), y=DataList[[i]][,1], cost=rtmspe,
			K=k, R=r, costArgs=list(trim=0.1))

	cvResults<-cvReshape(noIntCV)$reps
	cvResults<-cbind("noInt_CV"=cvResults[,"CV"], "Int_CV"=cvReshape(IntCV)$reps[,"CV"])
	
	cvR[[i]]<-cvResults
	}
return(cvR)

}
CV(TestCh1)
CV(TestCh2)
CV(TestPlasmid)


###########################################
# Function to plot proportions of gc or   #
# coding content for a specified binwidth #
# *Binwidths should be in data already    #
# Takes in the list of binwidths props    #
# Explorts graphs to current directory    #
# Might want to specify separate folders  #
# for different chromosomes               #
###########################################
PropVarPlot<-function(DataList, directory){
	for(i in 1:length(DataList)){
		InsertionFreq<-DataList[[i]][,1]
		gcProp<-DataList[[i]][,2]
		codingProp<-DataList[[i]][,3]

png(paste0("codingProp_", names(DataList)[i], ".png"))
		barplot(codingProp, xlab="Bins", ylab="Coding Region Percentage", cex.lab=1.4, xaxt='n')
dev.off()

png(paste0("gcProp_", names(DataList)[i], ".png"))
		barplot(gcProp[1:50], xlab="Bins", ylab="Percentage of GC-content", cex.lab=1.4, xaxt='n')
dev.off()


	#barplot(InsertionFreq, xlab="Bins", ylab="Number of insertions")
	}#End of for loop
}
