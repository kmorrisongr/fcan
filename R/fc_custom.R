# Collection of functions that are specific to your dataset

#' ED608 Convert Feature String To Nicer String
#'
#' A function that takes in an ugly Fc Array compound title and makes it nicer. Useful if you want axis labels/plot titles that aren't ugly.
#'
#' @param compound_feat A feature string that is ugly and contracted with excess periods, etc.
#'
#' @return The new feature string.
#'
#' @export
cust_ed608_pretty_title = function(compound_feat){
	reag = strsplit(compound_feat,"_")[[1]][1] 
	ant = strsplit(compound_feat,"_")[[1]][2] 
	feature_title = paste(gsub("[.]", "-", reag), gsub("[.]"," ",ant), sep="_")

	if (ant == "gp41.Ectodomain..HIV.1."){
		feature_title = paste(gsub("[.]", "-",reag), "gp41 Ectodomain HIV-1", sep="_")
	}

	if (ant == "P1.PE.synthetic.lipopeptide"){
		feature_title = paste(gsub("[.]", "-", reag), "P1-PE synth. lipopep.", sep="_")
	}

	return(feature_title)
}

# TODO: Update for subject column change

#' ED608 Bin Subjects By Challenge Information
#'
#' Regroup subjects based on whether the survived/remained uninfected for greater or less than the median number of challenges.
#'
#' @param fc The Fc Array data frame.
#' @param surv The survival information data frame. Should contain subject id, challenge to infection, and infection status.
#' @param new_groups The new group labels you wish to use.
#'
#' @return The Fc Array data frame with the group information replaced with the new group assignments.
#'
#' @export
cust_ed608_bin = function(fc, surv, new_groups=c("SUSC", "RES")){
	times = unique(fc$times)

	# At the first time point, all subjects should be present
	bins = fcff_only_attr(fc, times[1], times, "times")
	# Trim off the time point "A" at the end of the subject name
	rownames(bins) = substr(rownames(bins),1,nchar(rownames(bins))-1)
	# Get the survival data for the current data - in case we get passed a subset of fc
	local_surv = fccu_match_rows(surv$id,rownames(bins),surv)

	# Split at the median
	bins[(local_surv$challenges < median(local_surv$challenges)),1] = new_groups[1]
	bins[(local_surv$challenges >= median(local_surv$challenges)),1] = new_groups[2]

	# Assign the new groups to the whole dataset
	new_groups = vector("character",length=nrow(fc))
	for (i in 1:length(new_groups)){
		name = fc$subject[i]
		j = which(rownames(bins) == name)
		new_groups[i] = bins[j,1]
	}
	new_groups = as.factor(new_groups)
	new_groups = relevel(new_groups,ref=new_groups[1])

	group_clm = which(colnames(fc) == "group")
	fc$group = new_groups
	fc$group = as.factor(fc$group) # just in case
	colnames(fc)[group_clm] = "group"

	return(fc)
}

#' RV144 Bin Subjects By Challenge Information
#'
#' Regroup subjects based on whether the survived/remained uninfected for greater or less than the median number of challenges.
#'
#' @param fc The Fc Array data frame.
#' @param samples The sample mapping informatin data frame. Used for converting sample aliases (4.6.6, etc.) to PIDs, which is what survival information is paired with.
#' @param surv The survival information data frame. Should contain subject id, challenge to infection, and infection status.
#'
#' @return The Fc Array data frame with the group information replaced with the new group assignments.
#'
#' @export
cust_rv144_bin = function(fc, samples, surv){
	bins = vector("character")

	# Pull the PID from samples, use that to find infect status in surv
	for (s in fc$subject){
		pid = samples$PID[which(samples[,1] == s)]
		bins = c(bins, tolower(surv$infected[which(surv$id == pid)]))
	}

	fc$group = bins
	fc$group = as.factor(fc$group)

	return(fc)
}
