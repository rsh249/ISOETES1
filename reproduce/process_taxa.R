
##Get supplementary files from Pedersen et al. 2016
##
##Lake summary data are directly from:
#https://www.nature.com/nature/journal/v537/n7618/source_data/nature19085-sf6.xlsx
#and https://www.nature.com/nature/journal/v537/n7618/source_data/nature19085-sf5.xlsx
charlie = read.table('./reproduce/charlie_lake.csv', sep=',', header=TRUE)
silver = read.table('./reproduce/silver_lake.csv', sep=',', header=TRUE) #by silver lake I obviously mean spring lake

charlie=charlie[,1:58]
silver=silver[,1:77]
charlie=charlie[1:14,]
for (i in 2:ncol(charlie)){
	charlie[,i] = as.numeric(charlie[,i]);
}

ch.plant.reads = apply(charlie[,3:ncol(charlie)], 1, sum, na.rm=TRUE)

ch.reads.total = sum(apply(charlie[,3:ncol(charlie)], 1, sum, na.rm=TRUE))

silver=silver[1:18,]
for (i in 2:ncol(silver)){
	silver[,i] = as.numeric(silver[,i]);
}

si.plant.reads = apply(silver[,4:ncol(silver)], 1, sum, na.rm=TRUE)

si.reads.total = sum(apply(silver[,4:ncol(silver)], 1, sum, na.rm=TRUE))


ch.presence = charlie
ch.presence[,3:ncol(charlie)] = (charlie[,3:ncol(charlie)]>0)*1L
ch.pr = ch.presence[,c('sample_name', 'kyr.BP', 'Artemisia', 'Picea', 'Populus', 'Potamogeton', 'Salix', 
	'Aegilops', 'Hordeum', 'Helianthus', 'Fragaria', 
	'Brassica', 'Alfalfa', 'Fritillaria','Brachypodium')]

si.presence = silver
si.presence[,4:ncol(silver)] = (silver[,4:ncol(silver)]>0)*1L
si.pr = si.presence[,c('sample_name', 'Modelled.date..cal..BP.', 'Artemisia', 'Picea', 'Pinus',
	'Alnus', 'Salix',
	'Populus', 'Shepherdia', 
	'Elaeagnus', 'Brachypodium',
	'Agrostis', 'Aegilops', 'Elymus',
	'Hordeum', 'Leymus', 'Secale', 'Carex',
	'Artemisia', 'Rosa', 'Astragalus', 'Typha',
	'Ceratophyllum', 'Myriophyllum', 'Potamogeton')]
	
si.matrix = as.matrix(si.pr[,3:ncol(si.pr)])
rownames(si.matrix) = si.pr[,2]
heatmap(si.matrix, Rowv=NA, 
	Colv=NA, scale='none', 
	col =c('lightgrey', 'blue4'),
	main = 'Holi -- Spring Lake (common)'
)	
ch.matrix = as.matrix(ch.pr[,3:ncol(ch.pr)])
rownames(ch.matrix) = ch.pr[,2]
heatmap(ch.matrix, Rowv=NA, 
	Colv=NA, scale='none', 
	col =c('lightgrey', 'blue4'),
	main = 'Holi -- Charlie Lake (common)'
)	
#dev.off()

soil.dirs = list.files('soils.res')

#get time files:
kgen = matrix(ncol=200, nrow=length(soil.dirs));
kgen = as.data.frame(kgen)
cgen = matrix(ncol=200, nrow=length(soil.dirs));
cgen = as.data.frame(cgen)
bgen = matrix(ncol=500, nrow=length(soil.dirs));
bgen = as.data.frame(bgen);


for (i in 1:length(soil.dirs)){
	kgenera = paste('soils.res', soil.dirs[i], 'genera_th.csv', sep = '/')
	cgenera = paste('soils.res', soil.dirs[i], 'centgenera_th.csv', sep = '/')
	bgenera = paste('soils.res', soil.dirs[i], 'blastgen_th.csv', sep = '/')
	reads = paste('soils.res', soil.dirs[i], 'krakhits.o', sep = '/')
	if(file.exists(kgenera)){
		read.count = read.table(reads);
		krak = read.table(kgenera)
#		krak = subset(krak, krak[,2]>(read.count*0.01))
		kgen[i,1:(nrow(krak)+1)]=c(soil.dirs[i], as.character(krak$taxon))
		cent = read.table(cgenera)
		cgen[i,1:(nrow(cent)+1)]=c(soil.dirs[i], as.character(cent$taxon))
		blast = read.table(bgenera)
		bgen[i,1:(nrow(blast)+1)]=c(soil.dirs[i], as.character(blast[,1]))
	}
}

kgen[9,1] = 'ERR1560020';
cgen[9,1] = 'ERR1560020';
bgen[9,1] = 'ERR1560020';

library(reshape2)
kgen.melt <- melt(kgen, id.var='V1')
kgen.cast = dcast(kgen.melt, V1 ~ value)
cgen.melt <- melt(cgen, id.var='V1')
cgen.cast = dcast(cgen.melt, V1 ~ value)
bgen.melt <- melt(bgen, id.var='V1')
bgen.cast = dcast(bgen.melt, V1 ~ value)

##Merge with format from above for Charlie Lake and Spring Lake series. Make Heatmaps for common (expected) noncontaminant genera
kch = kgen.cast[which(kgen.cast[,1] %in% soil.dirs[1:14]),]
kch.pr = cbind(ch.presence[,1:2],kch[,c('V1', intersect(colnames(kch), colnames(ch.presence)))])
cch = cgen.cast[which(cgen.cast[,1] %in% soil.dirs[1:14]),]
cch.pr = cbind(ch.presence[,1:2],cch[,c('V1', intersect(colnames(cch), colnames(ch.presence)))])
bch = bgen.cast[which(bgen.cast[,1] %in% soil.dirs[1:14]),]
bch.pr = cbind(ch.presence[,1:2],bch[,c('V1', intersect(colnames(bch), colnames(ch.presence)))])

ksi = kgen.cast[which(kgen.cast[,1] %in% soil.dirs[15:32]),]
ksi.pr = cbind(si.presence[,1:2],ksi[,c('V1', intersect(colnames(ksi), colnames(si.presence)))])
csi = cgen.cast[which(cgen.cast[,1] %in% soil.dirs[15:32]),]
csi.pr = cbind(si.presence[,1:2],csi[,c('V1', intersect(colnames(csi), colnames(si.presence)))])
bsi = bgen.cast[which(bgen.cast[,1] %in% soil.dirs[15:32]),]
bsi.pr = cbind(si.presence[,1:2],bsi[,c('V1', intersect(colnames(bsi), colnames(si.presence)))])

kraknum = 199 - apply(is.na(kgen), 1, sum)
centnum = 199 - apply(is.na(cgen), 1, sum)
blastnum = 499 - apply(is.na(bgen), 1, sum)
gen.classif = cbind(soil.dirs, kraknum, centnum, blastnum)

gen.class.si = cbind(bsi.pr[,2], gen.classif[15:32,])
gen.class.ch = cbind(bch.pr[,2], gen.classif[1:14,])

gen.bp.si = rbind(cbind(gen.class.si[,1:3], rep("Kraken", 18)), 
	cbind(gen.class.si[,c(1,2,4)], rep("Centrifuge", 18)), 
	cbind(gen.class.si[,c(1,2,5)], rep("MegaBLAST", 18))
	)
gen.bp.si = as.data.frame(gen.bp.si)
colnames(gen.bp.si) = c('Age', 'Sample_ID', 'Num.Genera', 'Method')
gen.bp.si[,1] = factor(gen.bp.si[,1], levels=sort(unique(as.numeric(as.character(gen.bp.si[,1])))))
gen.bp.si[,3] = as.numeric(as.character(gen.bp.si[,3]))

gen.bp.ch = rbind(cbind(gen.class.ch[,1:3], rep("Kraken", 14)), 
	cbind(gen.class.ch[,c(1,2,4)], rep("Centrifuge", 14)), 
	cbind(gen.class.ch[,c(1,2,5)], rep("MegaBLAST", 14))
	)
gen.bp.ch = as.data.frame(gen.bp.ch)
colnames(gen.bp.ch) = c('Age', 'Sample_ID', 'Num.Genera', 'Method')
gen.bp.ch[,1] = factor(gen.bp.ch[,1], levels=sort(unique(as.numeric(as.character(gen.bp.ch[,1])))))
gen.bp.ch[,3] = as.numeric(as.character(gen.bp.ch[,3]))

hold.ch = cbind(gen.bp.ch, rep('Charlie Lake Series', 14))
hold.si = cbind(gen.bp.si, rep('Spring Lake Series', 18))
colnames(hold.ch) = c("Age", "Sample_ID", "Num.Genera", "Method", "Series")
colnames(hold.si) = c("Age", "Sample_ID", "Num.Genera", "Method", "Series")
gen.bp.all = rbind(hold.ch, hold.si)
gen.bp.all[,1] = factor(gen.bp.all[,1], levels=sort(unique(as.numeric(as.character(gen.bp.all[,1])))))

library(lattice)
library(RColorBrewer)

myColours <- brewer.pal(5,"Set2")

## Create your own list with
my.settings <- list(
  superpose.polygon=list(col=myColours[1:3], border="transparent"),
  strip.background=list(col=myColours[5]),
  strip.border=list(col="black")
)

tiff('soil_meta_gencount.tif', height =8, width=8, units='in', res=450);
lattice::barchart(Age~Num.Genera | Series,data=gen.bp.all, groups=Method,
         scales=list(x=list(rot=90,cex=0.8)), xlim=c(0,150), ylab = "Age (Calibrated Years BP)", main="Genera Classified in sedimentary aDNA series", 
         auto.key=list(space="top", columns=3, 
                       title="Method", cex.title=0.8), 
         par.settings = my.settings
)

dev.off()

##Heatmap comparison
##NOTES:
#differences in charlie lake ALL taxa
#differences in silver lake Populus, Picea, 


#NO Differences in silver lake Ceratophyllum in any method!

kch.matrix = as.matrix(kch.pr[,4:ncol(kch.pr)])
rownames(kch.matrix) = kch.pr[,2]
heatmap(kch.matrix, Rowv=NA,
        Colv=NA, scale='none',
        col =c('lightgrey', 'blue4'),
        main = 'Kraken -- Charlie Lake'
)
cch.matrix = as.matrix(cch.pr[,4:ncol(cch.pr)])
rownames(cch.matrix) = cch.pr[,2]
heatmap(cch.matrix, Rowv=NA,
        Colv=NA, scale='none',
        col =c('lightgrey', 'blue4'),
        main = 'Centrifuge -- Charlie Lake'
)

bch.matrix = as.matrix(bch.pr[,4:ncol(bch.pr)])
rownames(bch.matrix) = bch.pr[,2]
heatmap(bch.matrix, Rowv=NA,
        Colv=NA, scale='none',
        col =c('lightgrey', 'blue4'),
        main = 'MegaBLAST -- Charlie Lake'
)

ksi.matrix = as.matrix(ksi.pr[,4:ncol(ksi.pr)])
rownames(ksi.matrix) = ksi.pr[,2]
heatmap(ksi.matrix, Rowv=NA,
        Colv=NA, scale='none',
        col =c('lightgrey', 'blue4'),
        main = 'Kraken -- Spring Lake'
)
csi.matrix = as.matrix(csi.pr[,4:ncol(csi.pr)])
rownames(csi.matrix) = csi.pr[,2]
heatmap(csi.matrix, Rowv=NA,
        Colv=NA, scale='none',
        col =c('lightgrey', 'blue4'),
        main = 'Centrifuge -- Spring Lake'
)

bsi.matrix = as.matrix(bsi.pr[,4:ncol(bsi.pr)])
rownames(bsi.matrix) = bsi.pr[,2]
heatmap(bsi.matrix, Rowv=NA,
        Colv=NA, scale='none',
        col =c('lightgrey', 'blue4'),
        main = 'MegaBLAST -- Spring Lake'
)

##Function taxon/Lake specific heatmap
sitaxheat <- function(tax) {
	tax.mat = cbind(ksi.matrix[, tax], csi.matrix[,tax], bsi.matrix[,tax], si.matrix[,tax])
	rownames(tax.mat) = bsi.pr[,2];
	colnames(tax.mat) = c('Kraken', 'Centrifuge', 'MegaBLAST', 'Holi');
	par(cex.main=0.8) 
	heatmap(tax.mat, Rowv=NA,
      		Colv=NA, scale='none',
       		col =c('lightgrey', 'blue4'),
       		main = tax, cexRow=0.8, cexCol=0.7
	)
}


sitaxheat('Populus');
sitaxheat('Picea');
sitaxheat('Salix');
sitaxheat('Ceratophyllum');
sitaxheat('Potamogeton');

chtaxheat <- function(tax) {
	tax.mat = cbind(kch.matrix[, tax], cch.matrix[,tax], bch.matrix[,tax], ch.matrix[,tax])
	rownames(tax.mat) = bch.pr[,2];
	colnames(tax.mat) = c('Kraken', 'Centrifuge', 'MegaBLAST', 'Holi');
	par(cex.main=0.8) 
	heatmap(tax.mat, Rowv=NA,
      		Colv=NA, scale='none',
       		col =c('lightgrey', 'blue4'),
       		main = tax, cexRow=0.8, cexCol=0.7
	)
}

                                                                                                                           
chtaxheat('Populus');
chtaxheat('Picea');
chtaxheat('Salix');
chtaxheat('Potamogeton');

#chtaxheat('Ceratophyllum');
dev.off()

fulltaxheat <- function(tax) {
  tax.mat = matrix(data = 0, nrow=14, ncol=4);
  	#tax.mat = cbind(kch.matrix[, tax], cch.matrix[,tax], bch.matrix[,tax], ch.matrix[,tax])
    if(tax %in% colnames(kch.matrix)){
      tax.mat[,1] = kch.matrix[,tax]
    }
  if(tax %in% colnames(cch.matrix)){
    tax.mat[,2] = cch.matrix[,tax]
  }
  if(tax %in% colnames(bch.matrix)){
    tax.mat[,3] = bch.matrix[,tax]
  }
  if(tax %in% colnames(ch.matrix)){
    tax.mat[,4] = ch.matrix[,tax]
  }
  
	rownames(tax.mat) = bch.pr[,2];
	colnames(tax.mat) = c('Kraken', 'Centrifuge', 'MegaBLAST', 'Holi');
	
	tax.mat2 = matrix(data=0, nrow=18, ncol=4);
	#tax.mat = cbind(kch.matrix[, tax], cch.matrix[,tax], bch.matrix[,tax], ch.matrix[,tax])
	if(tax %in% colnames(ksi.matrix)){
	  tax.mat2[,1] = ksi.matrix[,tax]
	}
	if(tax %in% colnames(csi.matrix)){
	  tax.mat2[,2] = csi.matrix[,tax]
	}
	if(tax %in% colnames(bsi.matrix)){
	  tax.mat2[,3] = bsi.matrix[,tax]
	}
	if(tax %in% colnames(si.matrix)){
	  tax.mat2[,4] = si.matrix[,tax]
	}
	rownames(tax.mat2) = bsi.pr[,2];
	colnames(tax.mat2) = c('Kraken', 'Centrifuge', 'MegaBLAST', 'Holi');
	tax.mat = rbind(tax.mat, tax.mat2);
	tax.mat = tax.mat[order(as.numeric(rownames(tax.mat))),]
	return(tax.mat)
	par(cex.main=0.8) 
  heatmap(tax.mat, Rowv=NA,
      		Colv=NA, scale='none',
       		col =c('lightgrey', 'blue4'),
       		main = tax, cexRow=0.8, cexCol=0.7
	)

}


# tiff('soil_meta_picea.tif', height = 5, width = 4, units = 'in', res=450)
# fulltaxheat('Picea');
# dev.off()
# tiff('soil_meta_populus.tif', height = 5, width = 4, units = 'in', res=450)
# fulltaxheat('Populus');
# dev.off()
# tiff('soil_meta_potamogeton.tif', height = 5, width = 4, units = 'in', res=450)
# fulltaxheat('Potamogeton');
# dev.off()
# tiff('soil_meta_salix.tif', height = 5, width = 4, units = 'in', res=450)
# fulltaxheat('Salix');
# dev.off()
# tiff('soil_meta_hordeum.tif', height = 5, width = 4, units = 'in', res=450)
# fulltaxheat('Hordeum');
# dev.off()
# tiff('soil_meta_Ceratophyllum.tif', height = 5, width = 4, units = 'in', res=450)
# sitaxheat('Ceratophyllum');
# dev.off();

##Complex heatmaps:
library(pheatmap)
library(dplyr)


pinus = fulltaxheat("Pinus")
colnames(pinus) = paste(colnames(pinus), "Pinus", sep = '.')
picea = fulltaxheat("Picea")
colnames(picea) = paste(colnames(picea), "Picea", sep = '.')
populus = fulltaxheat("Populus")
colnames(populus) = paste(colnames(populus), "Populus", sep = '.')
Salix = fulltaxheat("Salix")
colnames(Salix) = paste(colnames(Salix), "Salix", sep = '.')
Ceratophyllum = fulltaxheat("Ceratophyllum")
colnames(Ceratophyllum) = paste(colnames(Ceratophyllum), "Ceratophyllum", sep = '.')
Potamogeton = fulltaxheat("Potamogeton")
colnames(Potamogeton) = paste(colnames(Potamogeton), "Potamogeton", sep = '.')
Artemisia = fulltaxheat("Artemisia")
colnames(Artemisia) = paste(colnames(Artemisia), "Artemisia", sep = '.')
Hordeum = fulltaxheat("Hordeum")
colnames(Hordeum) = paste(colnames(Hordeum), "Hordeum", sep = '.')
Elymus = fulltaxheat("Elymus")
colnames(Elymus) = paste(colnames(Elymus), "Elymus", sep = '.')


merged = cbind(populus, Salix, Artemisia, pinus, picea, Elymus, Hordeum, Potamogeton, Ceratophyllum)
methods = strsplit(colnames(merged), split="\\.")
for(i in 1:length(methods)){methods[[i]] = methods[[i]][1]}
methods=unlist(methods)
merged.factorDS = cbind(colnames(merged), 
                        c(rep("Populus", 4), rep("Salix", 4), rep("Artemisia", 4),
                          rep("Pinus", 4), rep("Picea", 4),
                          rep("Elymus",4), rep("Hordeum",4),
                          rep("Potamogeton",4), rep("Ceratophyllum", 4)), 
                        methods
                        )
colnames(merged.factorDS) = c('samplename', 'Genus', 'Method')


library(RColorBrewer)

myColours <- c(brewer.pal(9,"Paired"), brewer.pal(8,"Set1"))


genCols = myColours[1:9]
methCols = c("darkorange", "forestgreen", "blue2", "grey80")
names(genCols) = unique(merged.factorDS[,2])
names(methCols) = unique(merged.factorDS[,3])
AnnColour = list(Genus=genCols, Method=methCols)
merged.factorDS2 = merged.factorDS[,2:3]
rownames(merged.factorDS2) = merged.factorDS[,1]
merged.factorDS = as.data.frame(merged.factorDS2)
pheatmap(merged, color=rev(c('darkblue', 'lightgrey')),
         annotation_col=merged.factorDS, 
         annotation_colors=AnnColour, 
         cluster_cols=FALSE, 
         cluster_rows=FALSE, border_color = "black", gaps_col = c(4, 8, 12, 16, 20, 24,28, 32), 
         legend=TRUE, legend_breaks=c(0,1), legend_labels = c("Absent", "Present"),
         main='Ice Free Corridor aDNA classification of important taxa',
         ylab = "Years Before Present", show_colnames=FALSE, filename='ifc_heatmap.tiff'
         
         )



