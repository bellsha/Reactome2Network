#the purpose of this is to move everything into a format that allows the ids to easily go between species

desc<-unique(read.table("/Reactome/Input53/ReactomePathways.txt", sep='\t', comment.char='',na.strings='', quote="\""))
colnames(desc)<-c("ReactomeID","Name","Species")
#remove chicken...bc its difficult
desc<-subset(desc, Species != "Gallus gallus")
test<-reshape(desc, direction="wide", v.names="ReactomeID", idvar="Name", timevar="Species")
#note that warnings are thrown that signal that there are multiple IDs for a single mapping
mapedIds<-na.omit(c(as.matrix(test[,2:ncol(test)])))
#these occur bc there are places where the same "name" occurs in different branches
#but these have unique ids. but since the proteins are the same it *shouldnt* matter
unmap<-subset(desc, ReactomeID %in% setdiff(desc$ReactomeID,mapedIds))
#trying to reshape again this time without the ids that werent mapped
test2<-reshape(unmap, direction="wide", v.names="ReactomeID", idvar="Name", timevar="Species")
#need to get the columns the same
noname<-setdiff(colnames(test), colnames(test2))
Y<-as.data.frame(matrix(nrow=1,ncol=length(noname), NA)); colnames(Y)<-paste(noname)
x<-cbind(test2, Y)
name2ID<-rbind(test,x)
#save output
write.table(name2ID,"../Reactome/Output53/ReactomePathwayName2IDxSpec.txt", sep='\t', col.names=TRUE, quote=FALSE, row.names=FALSE)
