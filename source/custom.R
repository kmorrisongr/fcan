# Collection of functions that are specific to your dataset

source("fc_utils.R")

# A function for grabbing group IDs from subject number
cust_ruth_get_group = function(animal_num){
	if (animal_num <= 12) {
		return("M")
	} else if (animal_num <= 24){
		return("K")
	} else if (animal_num <= 36){
		return("L")
	}
}

# A function that extracts an animal's number from its string ID
cust_ruth_get_num = function(animal_id){
	# "TBRI 10A"/"TBRI 1A" -> "10"/" 1"
	id = substr(animal_id,nchar(animal_id)-2,nchar(animal_id)-1)
	# "10"/" 1" -> "10"/"","1"
	id = strsplit(id,' ')[[1]]
	# "10"/"","1" -> "10"/"1",""
	id = rev(id)
	# "10"/"1"
	return(as.integer(id[1]))
}

# A function for getting the color of a subject
cust_ruth_get_subj_col = function(group){
	if (group == "M"){
		return("black")
	} else if (group == "K"){
		return("blue")
	} else if (group == "L"){
		return("red")
	}
}

# A function for getting the shape of a subject
cust_ruth_get_subj_shape = function(group){
	if (group == "M"){
		return(21)
	} else if (group == "K"){
		return(22)
	} else if (group == "L"){
		return(23)
	}
}

# A function for getting the bin color of a subject
cust_ruth_get_bin_subj_col = function(group){
	if (group == "RES"){
		return("blue")
	} else if (group == "SUSC"){
		return("red")
	}
}

# A function for getting the bin shape of a subject
cust_ruth_get_bin_subj_shape = function(group){
	if (group == "RES"){
		return(21)
	} else if (group == "SUSC"){
		return(22)
	}
}

# A function that takes in an ugly Fc Array compound title and makes it nice
cust_ruth_pretty_title = function(compound_feat){
	reag = strsplit(compound_feat,"_")[[1]][1] 
	ant = strsplit(compound_feat,"_")[[1]][2] 
	feature_title = paste(gsub("[.]","-",reag), gsub("[.]"," ",ant), sep="_")

	if (ant == "gp41.Ectodomain..HIV.1."){
		feature_title = paste(gsub("[.]","-",reag), "gp41 Ectodomain HIV-1", sep="_")
	}

	if (ant == "P1.PE.synthetic.lipopeptide"){
		feature_title = paste(gsub("[.]","-",reag), "P1-PE synth. lipopep.", sep="_")
	}

	return(feature_title)
}

# TODO: Update for subject column change
# Recategorize the data based on survival outcome, splitting at the median
cust_ruth_bin = function(fc,surv,times,new_groups=c("SUSC","RES")){
	# At the first time point, all subjects should be present
	bins = fcu_only_attr(fc, times[1], times, "times")
	# Trim off the time point "A" at the end of the subject name
	rownames(bins) = substr(rownames(bins),1,nchar(rownames(bins))-1)
	# Get the survival data for the current data - in case we get passed a subset of fc
	local_surv = fcu_match_rows(surv$id,rownames(bins),surv)

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

cust_rv144_bin = function(fc,samples,surv){
	bins = vector("character")

	# Pull the PID from samples, use that to find infect status in surv
	for (s in fc$subject){
		pid = samples$PID[which(samples[,1] == s)]
		bins = c(bins,tolower(surv$infected[which(surv$id == pid)]))
	}

	fc$group = bins
	fc$group = as.factor(fc$group)

	return(fc)
}
