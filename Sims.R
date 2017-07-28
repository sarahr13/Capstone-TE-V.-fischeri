
###############################################
# Creating a simulated genome with insertions #
# Simulated insertion locations from a runif  #
# with the max and min as the beginning & end #
# of the chromosome			      #
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
# Function that replicates the simulation a specified number of times 			    #
# element: comes from the bin_calculator function -- ONLY input the bin width you want!     #
# elem_length: length of the chromosome/plasmid in bp					    #
# binwidth: the specified bin width; should be equal to the binwidth of the actual provided #
# numreps: how many replicates								    #
# numActual: number of insertions in the chromosome/plasmid				    #
# plotTitle: name of the title when it's output as a png				    #
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



