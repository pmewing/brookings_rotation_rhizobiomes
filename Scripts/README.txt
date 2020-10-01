Wrangle files if desired first. All scripts are in the ./Scripts/Wrangle folder. Execute in the order indicated by filenames.

Then, run the Overview and Part 1 .Rmd. 

Then, run Part 2 .Rmd. Consider switching the error estimation method from 'boot' to 'default', for SEM, and reducing the number of permutations. These options are in the first cell. SEM errors will change somewhat - the boot errors are more robust, but take a LONG time.

For both Rmd notebooks, note some of the commentary might be outdated, imprecise, or just wrong. A lot of it is copy-pasted.

Scripts should run just fine using packrat-managed dependencies. See documentation for packrat on RStudio's website.
