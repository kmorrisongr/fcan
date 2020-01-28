# ---------------------------
# EXPERIMENT 20200120
#	For (1) PLAC/VACC and (2) YES/NO amongst VACC:
# 		A. Perform univariate analysis for filtered dataset
#		B. Plot p-values and q-values
#		C. Get boxplots for features that are significant
# ---------------------------

experiment_name = "experiment_20200120"
config_location = "../fc_config_rv144.R"

source("../fc_load_config.R")
source("../../source/fc_uni.R")
source("../../source/fc_boxplot.R")


# ---------------------------
# EXPERIMENT-SPECIFIC SETUP
# ---------------------------

keep_filter = c("IgG3", "gp120", "gp140", "V1.V2.Clade.B", "V1.V2.Clade.E", "FcgRIIa", "MN")
discard_filter = c("C1q")

flags = list()
gopts = list()

flags$do_differs = FALSE
# Make univariate results comparable
flags$do_standard = TRUE
adj_method="fdr"
gopts$ppar_mar = c(5.1,4.1,4.1,6)
gopts$reag_legend_inset = c(-0.22,0)
gopts$box_legend_inset = c(-0.40,0)
gopts$bpar_mar=c(5.1,4.1,4.1,5.0)

source("funcs/exp_20200120.R")

# ---------------------------
# EXPERIMENT 1 - PLAC/VACC
# ---------------------------

flags$do_only_group = FALSE
flags$do_bin = FALSE
flags$nicer_legend = FALSE
gopts$bplot_group_cols = c("white","blue")
gopts$jplot_group_cols = c("black","blue")
gopts$lpp_ylim = NULL

# ---------------------------
# Initialize
final = fc_load_config(config_location, experiment_name)
fc = final$fc
results_dir = final$results_dir

exp_20200120(fc, results_dir, adj_method, flags, gopts)


# ---------------------------
# EXPERIMENT 2 - YES/NO VACC
# ---------------------------

flags$do_only_group = TRUE
flags$do_bin = TRUE
flags$nicer_legend = TRUE
gopts$bplot_group_cols = c("blue","red")
gopts$jplot_group_cols = c("blue","red")
gopts$lpp_ylim = c(0,0.8)

# ---------------------------
# Initialize
final = fc_load_config(config_location, experiment_name, only_group="VACCINE", bin_method=cust_rv144_bin)
fc = final$fc
results_dir = final$results_dir

exp_20200120(fc, results_dir, adj_method, flags, gopts)
