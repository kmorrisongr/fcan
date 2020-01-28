# ---------------------------
# PREAMBLE
# ---------------------------

# This file shows and comments on how to use what's present in the fcan repository
# Starting with fc_utils.R

source("source/fc_utils.R")


# ---------------------------
# INPUT & DATA CLEANUP
# ---------------------------

# Make sure data is in standard format!
	# Subject names are the row names
	# Group IDs in column 1
		# If not applicable for your data, fill fc[,"group"] with "dummy"
	# Time points in column 2
		# If not applicable for your data, fill fc[,"times"] with "dummy"
	# Feature IDs are the column names
fc = read.csv("example_data.csv",stringsAsFactors=FALSE,header=TRUE)
# Make sure you check that the fc loads correctly. I've had issues with rownames being in column 1 before when manually editing spreadsheets to fit the format

# Don't want any missing values!
# first_feat_col returns the first column that is not "times", "group", or "subject"
fc = fcu_impute_fc(fc, fcu_first_feat_col(fc),method=median)

# Grab treatment groups
groups = unique(fc[,"group"])
# Grab time points
# This is "dummy" for example_data.csv
times = unique(fc[,"times"])


# ---------------------------
# FEATURE FILTERING
# ---------------------------

# Remove highly correlated features if you like!
fc_less_cor = fcu_remove_cor(fc,0.95)

# Only look at these features
keep_filter = c("CH505")
# Alternatively, don't look at these features
discard_filter = c("p27")
# strict means only keep features that contain everything in keep_filter,
# 	or only discard features that contain everything in discard_filter
k_behavior = "permissive"
d_behavior = "permissive"
#k_behavior = "strict"

# Where do you want to output all our stuff? 
stats_dir = "./stats/"
dir.create(stats_dir)

# Wrapper that performs all the keyword filtering magic - make sure to save new directory name
filtered = fcu_wrap_filter(fc,stats_dir,keep_filter,discard_filter,k_behavior,d_behavior)
fc = filtered$fc
stats_dir = filtered$stats_dir

# ---------------------------
# FILL ME IN WITH EXAMPLE!
# Filter features based on whether they differ significantly from baseline or not
base_ids = (fc[,"group"] == "PLACEBO")


# ---------------------------
# GROUP/TIMEPOINT EXCLUSION
# ---------------------------

# Only look at same-side animals
fc_same = fcu_exclude_attr(fc,"separate","group")

# Suppose we had time point data
fc_bogus = rbind(fc,fc,fc)
fc_bogus[,"times"] = c(rep('A',nrow(fc)),rep('B',nrow(fc)),rep('C',nrow(fc)))

# Only data for time point A
# If we had real times generated above, we would just pass it here instead of c=(...)
fc_bogus_a = fcu_only_attr(fc_bogus, 'A', times=c('A','B','C'), "times")

# Want to get rid of only one?
fc_bogus_ac = fcu_exclude_attr(fc_bogus, 'B', "times")


# ---------------------------
# ANALYSIS PREP
# ---------------------------

# Make a list of datasets representing pairwise combinations of groups
# Add more groups for illustration
fc_bogus = rbind(fc,fc)
fc_bogus[,"group"] = c(rep("same",20),rep("separate",19),rep("lame",20),rep("seapirate",19))
fcs = fcu_fcs_combs(fc_bogus,groups=c("same","separate","lame","seapirate"))
