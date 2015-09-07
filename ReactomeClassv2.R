#######################################3
#Name: ReactomeClassv2.R
#Author: Shannon Bell shannon.bell2@nih.gov
#Purpose: this code prepares reactome annotation files for doing enrichment calcuations
#and returns a dataframe with the probe id along with other descriptors
#required packages: 
#required files from (http://www.reactome.org/pages/download-data/  ...version number added):
# ReactomePathwaysRelation.txt; ReactomePathways.txt; UniProt2Reactome.txt
#updated to try and improve the heirachery issues of one path being annotated to many different classes
#NOTE this version doesnt do nested network stuff
#there is a new file that gives the uniprot mapping at each level of the pathway heirchy...this might be useful for enrichment


######################################
#specify directories
inputDir<-"../Reactome/Input53"
outputDir<-"../Reactome/Output53"
#Resources

#reactome to uniprot
prot<-unique(read.table(file.path(inputDir, "UniProt2Reactome.txt"), sep='\t', comment.char='',na.strings='', quote="\""))
colnames(prot)<-c("UniProtID","ReactomeID","URL","Name","EvidenceCode","Species")
#reactome pathway relationships
relate<-unique(read.table(file.path(inputDir,"ReactomePathwaysRelation.txt"), sep='\t', comment.char='',na.strings='', quote="\""))
colnames(relate)<-c("parent","child")
#reactome pathway descriptions
desc<-unique(read.table(file.path(inputDir, "ReactomePathways.txt"), sep='\t', comment.char='',na.strings='', quote="\""))
colnames(desc)<-c("ReactomeID","Name","Species")

#Note that if you only want to work with one species,
#you can subset that from the desc file
#because some of the species have very different trees consider limiting it to select species that are appropriate for xcomparisons
# #(ie no CHICKEN bc it really messes things up with extra pathways)
# #basically getting a list of pathways that are not unique to chicken
# temp<-setdiff(unique(subset(desc, Species=="Gallus gallus")$Name),unique(subset(desc, Species!="Gallus gallus")$Name)); 
# validPaths<-setdiff(unique(subset(desc, Species != "Gallus gallus")$Name), temp)
#only using human
validPaths<-unique(subset(desc, Species %in% c("Homo sapiens", "Rattus norvegicus"))$Name)
#but it is easier to do it after the merge
#Each species has its own reactomID for a pathway name
t1<-merge(relate, desc, by.x="parent", by.y="ReactomeID", all.x=TRUE)
colnames(t1)<-c('parent','child','Parent.Name','Parent.Species')
data<-merge(t1, desc, by.x="child", by.y="ReactomeID", all.x=TRUE)
colnames(data)<-c('child','parent','Parent.Name','Parent.Species','Child.Name','Child.Species')
# #each species has its own distinct reactome ID
# human<-subset(data, Parent.Species == "Homo sapiens")
#extracting the parent-child relationships (removing species specific IDs)
pc<-unique(data[,c(3,5)]); pc<-pc[order(pc$Parent.Name),]
write.table(pc,file.path(outputDir,"ReactomePathwaysPC.txt"), sep='\t', col.names=TRUE, quote=FALSE, row.names=FALSE)
########################
########################
#to set up annotation for the Child classes using the name (no species specificity)
#this puts the top most level Major and will have the minor being no deeper thatn 4 levels down
#first removing paths that mess things up
pc2<-subset(pc, Child.Name %in% validPaths & Parent.Name %in% validPaths)
Child<-unique(pc2$Child.Name)
Parent<-unique(pc2$Parent.Name)
#these are the main classes that have no higher level 
notChild<-setdiff(Parent, Child)
relateN<-cbind(pc2, Major=NA,Minor=NA)
#assigns the highest level (ie parent name) as the major and minor name to those nodes which
#are not a child process of another node
for(i in 1:length(notChild)){
  relateN$Major[relateN$Parent.Name == notChild[i]]<-as.character(notChild[i])
  relateN$Minor[relateN$Parent.Name == notChild[i]]<-as.character(notChild[i])
}
#itterate through each parent->child relationship, 
#assigning the top main class as the "Major" class to each child term
#and one level beneath as the "Minor"
#nl is the number of levels down the tree to get the minor groupings
#these are useful for catagorization in cytoscape...at n=1 its just the same as Major
# on v50 2 results in 163 minor groupings while n=3 results in 378
nl<-2
#this is a counter for how getting to nl
n<-1
#while (n <4){
while (nrow(subset(relateN, is.na(Major)==TRUE)) >=1){
  #this brings the highest level as the "Major" pathways
  t2<-unique(subset(relateN, Major != "NA")[,c("Child.Name", "Major")])
  #assigns the highest level to the child term
  for(i in 1:nrow(t2)){
    relateN$Major[relateN$Parent.Name == t2$Child.Name[i]]<-as.character(t2$Major[i])
  }
  #using assignments from the last itteration (so n-1)

  if (n < nl){
    relateN$Minor<-relateN$Parent.Name
    relateN$Minor[is.na(relateN$Major)]<-NA
  }
  #once we have gone the desired depth, just propagate the same Major/Minor catagories down
  if(n >= nl){
   # y<-setdiff(relateN$Child.Name, t3$Child.Name)
    for(j in 1:nrow(t3)){
      relateN$Minor[relateN$Parent.Name == t3$Child.Name[j]]<-as.character(t3$Minor[j])
    }
    relateN$Minor[is.na(relateN$Major)]<-NA
  }
  #for getting annotation for the minor label
  t3<-unique(subset(relateN, Major != "NA" & Child.Name %in% setdiff(relateN$Child.Name, t2$Child.Name)))
  #this increases the level up one
  n<-n+1
}
relateN<-relateN[,c("Child.Name","Parent.Name","Major","Minor")]
#
#add in the top level annotations
temp<-cbind(Child.Name=notChild, Parent.Name=notChild, Major=notChild, Minor=notChild)
relateN<-rbind(temp, relateN)

#removal of the multiple cross labeling wrt major/minor
#this gets all the paths that are at the second level but under another pathway as a child
lev2<-unique(relateN$Minor)#this is all the level 1 and level 2
truChild<-setdiff(unique(relateN$Child.Name), lev2)
#now just select out the true child processes
relateNchild<-subset(relateN, Child.Name %in% truChild)
#Now need to add BACk in the top level annotation for those first 2 levels...
#these are the cases where the major=Minor=Parent and they are all top level
relateNupper<-subset(relateN, Parent.Name %in% notChild & Minor %in% notChild)
relateN2<-rbind(relateNchild, relateNupper)
#getting the merged annotation file
data2<-merge(desc, relateN2, by.x="Name", by.y="Child.Name", all=TRUE)
#taking care of the odd ball catagories (ones that werent in human and rat and appear to have annotated in an unstandard way)
missParent<-subset(data2, is.na(Parent.Name))
haveParent<-subset(data2, !is.na(Parent.Name))
temp<-merge(missParent[,c(1,2,3)], data[,c(1,3)], by.x="ReactomeID", by.y="child", all.x=TRUE)
temp2<-subset(temp, is.na(Parent.Name))
temp2$Parent.Name<-temp2$Name; newParent<-rbind(subset(temp, !is.na(Parent.Name)), temp2)
#for simplicity im going to tag this up even though it is likely breaking the structure
newParent$Minor<-newParent$Parent.Name
newParent$Major<-newParent$Parent.Name
data2b<-rbind(haveParent, newParent)
#this file gives the name, ID, species, parent name and major/minor pathway assiciation
write.table(data2b,file.path(outputDir, "ReactomePathwaysChildRelation.txt"), sep='\t', col.names=TRUE, quote=FALSE, row.names=FALSE)

#now to get the gene assignment
data3<-merge(prot, data2b, by=c("ReactomeID", "Name", "Species"), all=TRUE)
write.table(data3,file.path(outputDir,"ReactomePathways2UniProt.txt"), sep='\t', col.names=TRUE, quote=FALSE, row.names=FALSE)
rat<-subset(data3, Species =="Rattus norvegicus")
write.table(rat,file.path(outputDir,"ReactomePathways2UniProtRAT.txt"), sep='\t', col.names=TRUE, quote=FALSE, row.names=FALSE)

