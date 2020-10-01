# Rhizosphere microbial communities explain positive effects of diverse crop rotations on corn and soybean performance

Updated Oct 1st, 2020

Use this repository reproduce the analysis in:

Benitez-Ponce, M.S. et al, *in review*, Rhizosphere microbial communities explain positive effects of diverse crop rotations on corn and soybean performance. 

First, either `git clone` or download the repository and extract it to a folder of your choosing.

To view the unpolished results of the analysis without running, simply open the .html files in the ./Results folder. 

The repository is an RStudio project. Dependencies (e.g. `lavaan`) *should* be pretty well managed with packrat, and are in any case noted in the scripts (and will be auto-installed by many of scripts, as well). Therefore, once you've downloaded the repository, open the brookings_rotation_microbiome.Rproj file in RStudio. Then open the analysis files as desired. These are housed in ./Scripts. 

## Analysis Steps:
1. Wrangle files if desired first. All scripts are in the ./Scripts/Wrangle folder. Execute in the order indicated by filenames.

2. Run the Overview and Part 1 .Rmd notebook. Knit it to a new .html file using the button in RStudio.

3. Then, run Part 2 .Rmd. Consider switching the error estimation method from 'boot' to 'default', for SEM, and reducing the number of permutations. These options are in the first cell. SEM errors will change somewhat - the boot errors are more robust, but take a LONG time.

4. To make figures and export tables, run the scripts in ".Scripts/Figures and Tables".

### Notes
Paths should be maximally portable, thanks to the `here` package. 

Most variables (e.g. paths) are listed in the top few lines/code blocks of the script/notebook.

For both Rmd notebooks, note some of the commentary might be outdated, imprecise, or just wrong. A lot of it is copy-pasted. 

Writing new files is currently flagged out. Most scripts have a variable, either "ww" or "write", that is set to FALSE. Switch this to TRUE to enable writing. 

Send me a message with questions!