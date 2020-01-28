Purpose
-------
FCAN is a repository where you can find R functions for manipulating Fc Array datasets in a common format.

Common Format & Usage
-----
![Common Fc Array Format](https://github.com/kmorrisongr/fcan/blob/master/format.png)

Functions assume the structure of your data matches this. So if, for example, your data is not longitudinal, fill "times" with "dummy". Currently working on this not being a required component.

### Contents of this repository
This is how I have my directories laid out. If you construct a similar setup, you won't need to do much tweaking of fc\_config, and ideally none for fc\_load\_config.

* __parent/__
	* __example.R__ A work-in-progress file that shows examples of function usage
	* __results/__ contains the output of all your experiments
	* __source/__ contains all .R files (fc\_utils, fc\_boxplot, fc\_load\_config, ...)
	* __studies/__ contains all data (presumably in further sub-directories) for the studies you are working on
	* __working/__ contains fc\_config files, symlink (ln -s) to fc\_load\_config
		* __experiments/__ contains .R files that are your day-to-day "experiments"
			* __funcs/__ if you want to repeat an experiment (likely), build a function for it and stuff it here

Installation
------------
Clone this repository, or download fc\_utils.R somehow, and then source("fc\_utils.R") in your R code.

To-do
----
 * More documentation
 	* Expanded example file
	* More comments for each function
 * Get rid of any bugs that appear (just use ![XGH](https://gist.github.com/banaslee/4147370))

If something doesn't work like it should or you have something you wish was included here, send me an email kyle [dot] morrison [dot] gr [at] dartmouth [dot] edu
