# Rhizosphere microbial communities explain positive effects of diverse crop rotations on corn and soybean performance

Updated Oct 1st, 2020

Use this repository reproduce the analysis in:

Benitez-Ponce, M.S. et al, *in review*, Rhizosphere microbial communities explain positive effects of diverse crop rotations on corn and soybean performance. 

## Usage
First, either `git clone` or download the repository and extract it to a folder of your choosing.

To view the unpolished results of the analysis without running, simply open the .html files in the ./Results folder. 

The repository is an [RStudio](https://rstudio.com/) project. Dependencies (e.g. [`lavaan`](https://lavaan.ugent.be/)) *should* be pretty well managed with [packrat](https://rstudio.github.io/packrat/), and are in any case noted in the scripts (and will be auto-installed by many of scripts, as well, though maybe using a newer version). Therefore, once you've downloaded the repository, open the brookings\_rotation\_microbiome.Rproj file in RStudio. Then open the analysis files as desired. These are housed in ./Scripts. 

## File Structure

- "./Data" has the raw and wrangled data. It also has a number of .Rdata files, which contain results from the analysis and are accessed by the figure- and table-generating scripts.
- "./Results" has all figures and (unformatted) tables. It also has Rmarkdown notebooks.
- "./Scripts" are all scripts for running the analysis. Within this are a few files and folders:
	- "./Figures and Tables" has scripts for producing the figures/tables in the publication
	- "./Wrangle" has data wrangling scripts, to convert raw data to useable formats.
	- "./pcwOrd_0.2.1.Rdata" contains a snapshot of [`ordinator`](https://github.com/pmewing/ordinator) (f.k.a. `pcwOrd`) for running log-ratio analysis. Loading this will load all functions into your environment.
	- "./Convenience Functions.R" has a number of helper functions for Part 1, and a number which may not be used.
	- "./Part 2 - Helper Functions.R" similarly has a number of helper functions for Part 2, plus some that might not be used.
	- "./*.Rmd" are the Rmarkdown notebooks in which the analyses were run. 
- "packrat" contains the R libraries/versions used. 


## Analysis Steps:
1. Open "./brookings\_rotation\_microbiome.Rproj" in Rstudio. This will set your home directory and give access to the dependencies.
2. Wrangle files if desired first. All scripts are in ./Scripts/Wrangle. Execute in the order indicated by filenames. See notes.
2. Run the Overview and Part 1 .Rmd notebook. Knit it to a new .html file using the button in RStudio.
3. Then, run Part 2 .Rmd. Consider switching the error estimation method from 'boot' to 'default', for SEM, and reducing the number of permutations. These options are in the first cell. SEM errors will change somewhat - the boot errors are more robust, but take a LONG time.
4. To make figures and export tables, run the scripts in ".Scripts/Figures and Tables".

### Notes
1. Paths should be portable, thanks to the `here` package.
2. Most variables (e.g. paths) are listed in the top few lines/code blocks of the script/notebook.
3. For both Rmd notebooks, note some of the commentary might be outdated, imprecise, or just wrong. A lot of it is copy-pasted. 
4. Writing new files is currently flagged out. Most scripts have a variable, either "ww" or "write", that is set to FALSE. Switch this to TRUE to enable writing. This variable is usually near the top, but might not be.
6. If running [FUNGuild](https://github.com/UMNFuN/FUNGuild), you'll need to dowload that an place Guilds_v1.1.py in a folder, "./FUNGuild". See the section `#### Update funguild` in "4 - Wrangle - ITS.R". If you're on windows you might need to change the script a bit to execute the external python command.
5. Send me a message with questions!
