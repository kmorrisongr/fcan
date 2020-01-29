# ---------------------------
# FILE LOCATIONS
# ---------------------------

# Where are all the .R files located?
source_dir = "../source/"

fc_location = "../studies/rv144/data_2014/fc_800_trimmed.csv"

# RV144-specific; used for converting sample aliases to PIDs (rv144_bin)
samples_location = "../studies/rv144/data_2014/sample_ids_800.csv"
samples = read.csv(file=samples_location,stringsAsFactors=FALSE,header=TRUE)

surv_location = "../studies/rv144/data_2014/survival_800.csv"

results_dir = "../results/"


# ---------------------------
# FEATURE FILTERING
# ---------------------------

if (!exists("keep_filter")){
	keep_filter = c()
}
if (!exists("discard_filter")){
	discard_filter = c()
}

k_behavior = "permissive"
d_behavior = "permissive"

if (!exists("flags")){
	flags = list()
}


# ---------------------------
# PARAMETERS
# ---------------------------

# Are experiments one directory level below fc_config.R?
flags$adj_cwd = TRUE

flags$do_differs = FALSE

# ---------------------------
# GRAPHICS OPTIONS
# ---------------------------

if (!exists("gopts")){
	gopts = list()
}
gopts$time_cols = c("black")

# For everything not in feat_ant, you'll get the last ant_shape
gopts$feat_ant = c("gp120","gp140","V1.V2.","gag","p17,p24,p51,p55,p66")
gopts$ant_shapes = c(21,22,23,24,25,7)
gopts$feat_ant_legend = c("gp120","gp140","V1.V2.","gag","p-antigens")

# Put the least specific (like IgG) last to break out by subclass (really, by substring)
	# Put it first if you want them collapsed into a single category
	# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# For everything not in these categories, you'll get the last reag_col
gopts$reag_cats = c("IgG1","IgG2","IgG3","IgG4","FcgRIIa","FcgRIIb","FcgRIIIa","FcgRIIIb","C1q","IgG")
gopts$reag_cols = c("blue","pink","turquoise1","saddlebrown","forestgreen","purple","gold","orange","green","black","orchid")

# How we decide what the categories for boxplots, longitudinal plots, etc. are
gopts$plot_antigens = c("gp41.HxBc2","gp41.Ectodomain..HIV.1.","P1.PE.synthetic.lipopeptide")

gopts$volc_xlim = c(-2,2)
gopts$volc_ylim = c(0,4)
gopts$reag_legend_inset = c(-0.20,0) 
