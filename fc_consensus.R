# Written by Kyle Morrison for analysis of Fc Array Data
# This function figures out consensus features based on the output of fc_foch_pvals()

fc_consensus = function(stats_dir,cons_gen_method,cons_cut){
	cons = read.table(file=paste(stats_dir,cons_gen_method,"_consensus_sig_feats.csv",sep=""),
			  sep=",",stringsAsFactors=FALSE,header=TRUE)

	# Get occurrences of all features, removes the NA tally automagically
	cons_occ = cons_plot = as.data.frame(table(unlist(cons)))
	
	# Write out those that are present in >= cons_cut of time point comparisons
	keep = which(cons_occ[,2] >= ncol(cons)*cons_cut)
	cons_occ = cons_occ[keep,]
	order = order(cons_occ[,2],decreasing=TRUE)
	cons_occ = cons_occ[order,]
	write.table(cons_occ,file=paste(stats_dir,"cons_final.csv",sep=""),
		    sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

	# Plotting frequency vs feature
	png(file=paste(stats_dir,"cons_plot.png",sep=""))
	plot(cons_plot[,2],xlab="Feature",ylab="Frequency",pch=21,col="black",bg="black")
	abline(h=ncol(cons)*cons_cut,col="red")
	dev.off()
}
