########################################################
#     Transform FASTA file to vector of 0's and 1's    #
# Can also calculate the AT and GC content proportions #
#      Just enter in fasta file name including ""      #
#         Note: FASTA file needs to be in pwd          #
########################################################
Transform<-function(fastaname){
	library(seqinr)
	Sequence<-read.fasta(fastaname)
	Seq<-unlist(Sequence)

  #Replacing w/ 0's and 1's
  Seq<-replace(Seq, Seq=="a", 0)
  Seq<-replace(Seq, Seq=="t", 0)
  Seq<-replace(Seq, Seq=="g", 1)
  Seq<-replace(Seq, Seq=="c", 1)

  AT<-length(which(Seq==0)) #Total AT content
  GC<-length(which(Seq==1)) #Total GC content

  ATprop<-AT/length(Seq) 
  GCprop<-GC/length(Seq)

  Results<-list("seq"=Seq,"ATprop"=ATprop, "GCprop"=GCprop)
  return(Results) #Only returns the AT and GC proportions 
}#End of function


#Performing on chromosome 1
Ch1vec<-Transform("Ch1.fasta")
Ch1vec$ATprop
Ch1vec$GCprop
write.table(Ch1vec$seq, file="Ch1_GCcontent.txt", row.names=F)

#Performing on chromosome 2
Ch2vec<-Transform("Ch2.fasta")
Ch2vec$ATprop
Ch2vec$GCprop
write.table(Ch2vec$seq, file="Ch2_GCcontent.txt", row.names=F)

#Plasmid
Pvector<-Transform("plasmid.fasta")
row.names(Pvector$seq)
write.table(Pvector$seq, file="Plasmid_GCcontent.txt", row.names=F)


#####################
# Different package #
#####################
source("https://bioconductor.org/biocLite.R") 
biocLite("Biostrings")
library(Biostrings)  
p<-readDNAStringSet("plasmid.fasta")
GC<-c("G", "C", "A", "T")
letterFrequency(p, GC)
