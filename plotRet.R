setwd("/beegfs/group_dv/home/RCui/killifish_genomes/miRNA");
datGenomeSize <- read.table("genomesize.txt", header = T, sep="\t");
datGenomeSize <- read.table("genomesize_more.txt", header = T, sep="\t");

nMismatch <- 2;
sFiles <- paste("out_maxdiff", nMismatch, "/*/count_for_*.txt", sep="");
arrF <- Sys.glob(sFiles);

datCounts <- NULL;
for(sF in arrF) {
  datC <- read.table(sF, header=T, sep="\t");
  datCounts <- rbind(datCounts, datC);
}

datM <- merge(datGenomeSize, datCounts, by.x="Code", by.y="Taxon", all = T);
datM <- datM[complete.cases(datM), ]
datM$Genus <- factor(datM$Genus, c('Kryptolebias', 'Austrofundulus', 'Aplocheilus', 'Pachypanchax', 'Epiplatys', 'Archiaphyosemion', 'Scriptaphyosemion', 'Callopanchax', 'Aphyosemion', 'Fundulopanchax', 'Fundulosoma', 'Pronothobranchius', 'Nothobranchius') );

datM$SeqCoverage <- datM$ProccessedBases / datM$Genomesize_BUSCO;
datM$PerHaploidCount <- datM$TotalHits / datM$SeqCoverage;

datSum <- aggregate(datM, by = list(Code=datM$Code) , FUN=function(x) {if (is.numeric(x)) {sum(x)} else {x[1]} } )
datSum$ProbeID <- 'Sum';
datSum$ProbeSeq <- 'All sequences summed';

datM <- rbind(datM, datSum[2:ncol(datSum)]);

write.table(datM, file=paste('mirna_search_mismatch_', nMismatch,"bp.txt", sep="") , col.names = T, row.names = F, quote=F, sep="\t");


pdf(file=paste('mirna_search_mismatch_', nMismatch,"bp.pdf", sep=""), width=8,height = 6 )
par(mar=c(6,8,4,1)+1)

arrSeq <- unique(datM$ProbeSeq);
for(sSeq in arrSeq) {
  datThisM <- datM[datM$ProbeSeq == sSeq, ];
  arrColors <- c('blue','red', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'blue', 'orange', 'orange', 'red', 'red' );
  boxplot(datThisM$PerHaploidCount ~ datThisM$Genus, main=paste(sSeq, "mismatch: ", nMismatch ) , col=arrColors, boxfill=arrColors, horizontal=TRUE,las=1, xlab="Counts per haploid genome" );
}

dev.off();