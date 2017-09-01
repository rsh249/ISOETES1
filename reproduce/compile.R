##R script to collect and summarize plant metagenomic results from Pedersen et al. 2016 data and simulations
##Run from base directory. (i.e., that one that has the soils.res folder in it)
## Pedersen et al. 2016 results

soil.dirs = list.files('soils.res')

#get time files:
t = matrix(ncol=3, nrow=length(soil.dirs));
t = as.data.frame(t)
counts = matrix(ncol=1, nrow=length(soil.dirs));
counts = as.data.frame(counts)
for (i in 1:length(soil.dirs)){
	timefile = paste('soils.res', soil.dirs[i], 'timefile.o', sep = '/')
	readcfile= paste('soils.res', soil.dirs[i], 'read.count', sep = "/")
	if(file.exists(timefile)){
		time=read.table(timefile)
		reads=read.table(readcfile)
		t[i,] = time[1,1:3];
		counts[i,] = reads[1,1];
	}
}

colnames(t) = c('centrifuge', 'kraken', 'blast');

tiff('soil_meta_time.tif', height=6, width = 6, units='in', res = 350)
plot(counts[,1], log(t[,1]*4), ylim = c(0, 11), ylab='log[time(s)]', xlab='Number of reads')
points(counts[,1], log(t[,2]*4), pch = 15)
points(counts[,1], log(t[,3]), pch =10)
legend(20000000, 2, c('Centrifuge', 'Kraken', 'MegaBLAST'), pch = c(1, 15, 10))
dev.off()



#reads classified vs. read count
class = matrix(ncol=3, nrow=length(soil.dirs));
class = as.data.frame(class);

for (i in 1:length(soil.dirs)){
	classfile = paste('soils.res', soil.dirs[i], 'percentclass.o', sep = '/')
	if (file.exists(classfile)){
		if (file.size(classfile) > 1) {
			cl=read.table(classfile);
			print(cl);
			class[i,1:3] = cl[1,1:3];
		}
	}
}
colnames(class) = c('Centrifuge', 'Kraken', 'MegaBLAST');
#tiff('soil_meta_classrate.tif', height=6, width = 6, units='in', res = 350)
#plot(counts[,1], class[,1]*100, ylim =c(0,17), ylab='% of reads classified', xlab='Number of reads')
#points(counts[,1], class[,2]*100, pch = 15)
#points(counts[,1], class[,3]*100, pch =10)
#legend(23000000, 16, c('Centrifuge', 'Kraken', 'MegaBLAST'), pch = c(1, 15, 10))
#dev.off()

#Do bar chart or box/whiskers to show mean % classification rate for Pedersen data
tiff('soil_meta_classrate.tif', height=6, width = 6, units='in', res = 350)
boxplot(class*100, ylab='% of reads classified', xlab='Number of reads')
legend(23000000, 16, c('Centrifuge', 'Kraken', 'MegaBLAST'), pch = c(1, 15, 10))
dev.off()



###Simulation results
sim.dirs = list.files('meta_sim')
gtpr = matrix(ncol=3, nrow=length(sim.dirs));
gtpr = as.data.frame(gtpr);
otpr = matrix(ncol=3, nrow=length(sim.dirs));
otpr = as.data.frame(otpr);
gppp = matrix(ncol=3, nrow=length(sim.dirs));
gppp = as.data.frame(gppp);
oppp = matrix(ncol=3, nrow=length(sim.dirs));
oppp = as.data.frame(oppp);

groups = c(rep('10', 5), rep('100', 5), rep('50', 5))
for (i in 1:length(sim.dirs)){
	gtprfile = paste('meta_sim', sim.dirs[i], 'gtpr.o', sep = '/')
        otprfile = paste('meta_sim', sim.dirs[i], 'otpr.o', sep = '/')
        gpppfile = paste('meta_sim', sim.dirs[i], 'gppp.o', sep = '/')
        opppfile = paste('meta_sim', sim.dirs[i], 'oppp.o', sep = '/')
	if(file.exists(opppfile)){
		gtpr[i,1:3] = read.table(gtprfile)[,1:3]
		otpr[i,1:3] = read.table(otprfile)[,1:3]
		gppp[i,1:3] = read.table(gpppfile)[,1:3]
		oppp[i,1:3] = read.table(opppfile)[,1:3]
	
	}

}

gtpr = cbind(groups, gtpr)
otpr = cbind(groups, otpr)
gppp = cbind(groups, gppp)
oppp = cbind(groups, oppp)
colnames(gtpr) = c('group', 'Kraken', 'Centrifuge', 'MegaBLAST')
colnames(otpr) = c('group', 'Kraken', 'Centrifuge', 'MegaBLAST')
colnames(gppp) = c('group', 'Kraken', 'Centrifuge', 'MegaBLAST')
colnames(oppp) = c('group', 'Kraken', 'Centrifuge', 'MegaBLAST')


gtpr.m=aggregate(gtpr[,2:4], by=list(gtpr$group), mean, na.rm=TRUE)
otpr.m=aggregate(otpr[,2:4], by=list(otpr$group), mean, na.rm=TRUE)
gppp.m=aggregate(gppp[,2:4], by=list(gppp$group), mean, na.rm=TRUE)
oppp.m=aggregate(oppp[,2:4], by=list(oppp$group), mean, na.rm=TRUE)


##NEXT PART DOES NOT WORK:
gtpr.m = gtpr.m[c(1,3,2),1:4]
gtpr.m[,1] = as.numeric(as.character(gtpr.m[,1]))


otpr.m = otpr.m[c(1,3,2),1:4]
otpr.m[,1] = as.numeric(as.character(otpr.m[,1]))



gppp.m = gppp.m[c(1,3,2),1:4]
gppp.m[,1] = as.numeric(as.character(gppp.m[,1]))


oppp.m = oppp.m[c(1,3,2),1:4]
oppp.m[,1] = as.numeric(as.character(oppp.m[,1]))


tiff('sim_classifiers.tiff', height = 12, width = 12, units='in', res =450)
par(mfrow=c(2,2))


plot(gtpr.m[,c(1,2)], type = 'p',axes=FALSE, pch =1, col=NA, ylim=c(0,1), ylab='Sensitivity', xlab = "ntax", main='Genus classification - True Positive Rate')
points(gtpr.m[,c(1,3)], type = 'b', pch = 15)
points(gtpr.m[,c(1,2)], type = 'b', pch =1)
points(gtpr.m[,c(1,4)], type = 'b', pch = 10)
axis(side = 1, at = c(10,50,100))
axis(side = 2, at = c(0,0.20, 0.40, 0.60, 0.80, 1))
legend(65, 0.40, c('Kraken', 'Centrifuge', 'MegaBLAST'), pch = c(1, 15, 10))


plot(otpr.m[,c(1,2)], type = 'p',axes=FALSE, pch =1, col=NA, ylim=c(0,1), ylab='Sensitivity', xlab = "ntax", main='Order classification - True Positive Rate')
points(otpr.m[,c(1,3)], type = 'b', pch = 15)
points(otpr.m[,c(1,2)], type = 'b', pch =1)
points(otpr.m[,c(1,4)], type = 'b', pch = 10)
axis(side = 1, at = c(10,50,100))
axis(side = 2, at = c(0,0.20, 0.40, 0.60, 0.80, 1))
legend(65, 0.40, c('Kraken', 'Centrifuge', 'MegaBLAST'), pch = c(1, 15, 10))


plot(gppp.m[,c(1,2)], type = 'p',axes=FALSE, pch =1, col=NA, ylim=c(0,1), ylab='Precision', xlab = "ntax", main='Genus classification - Positive Predictive Value')
points(gppp.m[,c(1,3)], type = 'b', pch = 15)
points(gppp.m[,c(1,2)], type = 'b', pch =1)
points(gppp.m[,c(1,4)], type = 'b', pch = 10)
axis(side = 1, at = c(10,50,100))
axis(side = 2, at = c(0,0.20, 0.40, 0.60, 0.80, 1))
legend(65, 0.40, c('Kraken', 'Centrifuge', 'MegaBLAST'), pch = c(1, 15, 10))



plot(oppp.m[,c(1,2)], type = 'p',axes=FALSE, pch =1, col=NA, ylim=c(0,1), ylab='Precision', xlab = "ntax", main='Order classification - Positive Predictive Value')
points(oppp.m[,c(1,3)], type = 'b', pch = 15)
points(oppp.m[,c(1,2)], type = 'b', pch =1)
points(oppp.m[,c(1,4)], type = 'b', pch = 10)
axis(side = 1, at = c(10,50,100))
axis(side = 2, at = c(0,0.20, 0.40, 0.60, 0.80, 1))
legend(65, 0.40, c('Kraken', 'Centrifuge', 'MegaBLAST'), pch = c(1, 15, 10))

dev.off()


