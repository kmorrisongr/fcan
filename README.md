Purpose
-------
FCAN is a repository where you can find R functions for manipulating Fc Array datasets in a common format.

Common Format & Usage
-----
![Common Fc Array Format](https://github.com/kmorrisongr/fcan/blob/master/format.png)

Functions assume the structure of your data matches this. Having longitudinal information ("times" column) is optional - if your data is single-time-point, no "times" column necessary.

### Recommended Layout
This is the suggested directory structure. Modify flags$adj_cwd if __experiments/__ does not live below fc\_load\_config.R

* __parent/__
	* __results/__ contains the output of all your experiments
	* __studies/__ contains all data (presumably in further sub-directories) for the studies you are working on
	* __working/__ contains fc\_config files
		* __experiments/__ contains .R files that are your day-to-day "experiments"
			* __funcs/__ if you want to repeat an experiment (likely), build a function for it and stuff it here

Installation
------------
In R:

install.packages("devtools")

library(devtools)

install\_github("kmorrisongr/fcan")

library(fcan)

Then you're off to the races!

To-do
----
 * More documentation
 	* example.R/wiki of some sort
 * Get rid of any bugs that appear (just use [XGH](https://gist.github.com/banaslee/4147370))

If something doesn't work like it should or you have something you wish was included here, send me an email kyle [dot] morrison [dot] gr [at] dartmouth [dot] edu
