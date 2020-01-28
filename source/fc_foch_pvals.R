# Written by Kyle Morrison for analysis of Fc Array Data
# This function performs statistical tests on the data between groups
	# Generates p-values and plots them
	# Generates significant feature tables

fc_foch_pvals = function(fcs, times, stats_dir, adj_method, line_method, flags, gopts, alternative="two.sided", start_clm=fcu_first_feat_col(fcs[[1]])){
	# ---------------------------
	# Do all the statistics!
	# ---------------------------

	for (i in 1:length(fcs)){
		fc_current = fcs[[i]]
		current_name = names(fcs)[i]

		fc_dir = paste(stats_dir,names(fcs)[i],"/",sep='')
		dir.create(fc_dir)
		
		foch_dir = paste(fc_dir,"fold_changes/",sep='')
		dir.create(foch_dir)

		master_ps = as.data.frame(matrix(nrow=ncol(fc_current[,-(1:(start_clm-1))]),ncol=length(times)))
		rownames(master_ps) = colnames(fc_current[,-(1:(start_clm-1))])
		colnames(master_ps) = times

		cat("working on",current_name,'\n')

		for (point in times){
			# Check for actual time points
			if (point != "dummy"){
				fc_point = fcu_only_attr(fc_current,point,times,"times")
				
			} else {
				fc_point = fc_current
			}
			
			p_vals = vector("numeric",length=( ncol(fc_point) - (start_clm-1) ))
			
			# Data frame for single comparison retains colname, for clarity when write()ing
			foch = as.data.frame(matrix(0,nrow=( ncol(fc_point) - (start_clm-1) ),ncol=1))
			colnames(foch) = current_name
			rownames(foch) = colnames(fc_point)[-(1:(start_clm-1))]

			# For each feature
			for (j in start_clm:ncol(fc_point)){
				local_data = data.frame(as.numeric(fc_point[,j]),fc_point[,1])
				colnames(local_data) = c("measure","groups")

				result = fcu_sing_comp_foch(local_data,current_name,
							method=median,alternative=alternative,exact=FALSE)

				foch[j-(start_clm-1),1] = result$sing_foch
				p_vals[j-(start_clm-1)] = result$p_val
			}

			if (adj_method != "none"){
				cat("adjusting p-values using",adj_method,'\n')
				p_vals = p.adjust(p_vals,method=adj_method)
			}

			master_ps[,point] = p_vals
			cat("number of features with p < 0.05:", length(which(p_vals < 0.05)), '\n')

			# Output fold change data
			file_name = paste(foch_dir, point, "_fold_change.csv", sep='')
			write.table(foch, file=file_name, sep=',', quote=FALSE, row.names=TRUE, col.names=TRUE)
		}

		# Write out all the p_vals to a file
		file_name = paste(fc_dir, adj_method, "_p_values.csv", sep='')
		write.table(master_ps, file=file_name, sep=',', quote=FALSE, row.names=TRUE, col.names=TRUE)


		# ---------------------------
		# Plot p-values

		png_name = paste(fc_dir,adj_method,"_p_values_plot.png",sep='')
		png(file=png_name)
		par(mar=c(5.1,4.1,4.1,6))
		plot(na.omit(master_ps[,1]), ylim=c(0,1),
		     xlab="Feature index", ylab="p-value", col=gopts$time_cols[1],
		     main=paste("fc", current_name, ":", adj_method, "p-values for each feature\n(by time point)", sep=' '))
		par(xpd=TRUE)
		legend("topright", inset=c(-0.15, 0), legend=times, bty="n", text.col=gopts$time_cols, col=gopts$time_cols, pt.bg=gopts$time_cols, 
		       cex=0.7)
		par(xpd=FALSE)

		if (ncol(master_ps) > 1){
			for (j in 2:ncol(master_ps)){
				points(na.omit(master_ps[,j]), col=gopts$time_cols[j])
			}
		}

		# If we adjust the values, we keep the default cutoff line
		if (adj_method != "none"){
			abline(h=0.05, col="red", lty=2)

		# Else if we don't adjust the values, maybe we move the cutoff line?
		} else {
			if (line_method == "bonferroni"){
				abline(h=(0.05/nrow(master_ps)), col="red", lty=2)

			} else if (line_method == "raw"){
				abline(h=0.05, col="red", lty=2)
			}
		}
		dev.off()


		# ---------------------------
		# Significant findings saving

		# Done this way in case you only have one column (time point)
		original_ps = master_ps
		significant = master_ps < 0.05
		keep = as.vector(apply(significant, 1, any, na.rm=TRUE))
		master_ps = as.data.frame(master_ps[keep,])
		rownames(master_ps) = rownames(original_ps)[keep]
		colnames(master_ps) = colnames(original_ps)
		
		if (nrow(master_ps) > 0){
			# What's the most p-values < 0.05 among all the time points?
			lengths = apply(master_ps, 2, function(x) { return( length(which(x < 0.05)) ) } )
			max_len = max(lengths)

			sig = as.data.frame(matrix(ncol=ncol(master_ps), nrow=max_len))
			colnames(sig) = times

			# For saving the ones we will modify then slam together
			sig_list = list()
			sig_names = vector("character")
			l = 1
			for (j in c("features","p_values")){
				if (j == "features"){
					# For each set of observations
					for (k in 1:ncol(master_ps)){
						to_write = rownames(master_ps)[which(master_ps[,k] < 0.05)]
						to_write = fcu_extend(to_write, max_len, NA)

						sig[,k] = to_write
					}

					sig_list[[l]] = sig
					sig_names = c(sig_names,names(fcs)[i])
					l = l + 1

				} else {
					for (k in 1:ncol(master_ps)){
						to_write = master_ps[which(master_ps[,k] < 0.05),k]
						to_write = fcu_extend(to_write,max_len,NA)

						sig[,k] = to_write
					}
				}

				file_name = paste(fc_dir, adj_method, "_significant_", j, ".csv", sep='')
				write.table(sig, file=file_name, sep=',', quote=FALSE, row.names=FALSE, col.names=TRUE)
			}
		}
	}


	# ---------------------------
	# Merge the feature tables
	# ---------------------------

	if (flags$do_cons){
		# Find the max number of rows
		max = 0
		for (i in 1:length(sig_list)){
			len = nrow(sig_list[[i]])	
			if (len > max){
				max = len
			}
		}

		new_sig_list = list()
		for (i in 1:length(sig_list)){
			to_add = sig_list[[i]]
			# So we can differentiate data in the mungo data.frame
			colnames(to_add) = paste( rep(sig_names[i],ncol(to_add)), colnames(to_add), sep='-')

			# If you are shorter than the longest, fix that
			if (nrow(to_add) < max){
				for (j in 1:(max-nrow(to_add))){
					to_add = rbind(to_add, NA)
				}
			}

			# If you have column/s made up of only NAs, nix that
			blank_count = as.vector(apply(to_add, 2, function(x){ length(which(is.na(x))) }))
			keep = which(blank_count != nrow(to_add))
			to_add = to_add[,keep]

			new_sig_list[[i]] = to_add
		}

		# Smash them all together, write it out
		mungo_sigs = as.data.frame(lapply(new_sig_list, cbind))
		write.table(mungo_sigs, file=paste(stats_dir, adj_method, "_consensus_sig_feats.csv", sep=''), 
			    sep=',', quote=FALSE, row.names=FALSE, col.names=TRUE)
	}
}
