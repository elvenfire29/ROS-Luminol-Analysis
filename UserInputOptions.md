# User Input Options

## User Input I: Steps 1-3 below must be changed to according to a specific dataset on a local computer.

Step 1:  Select the working directory. Indicate where the input file resides. The same folder will be the location for the output files.

Step 2: Select the input file. The input file should be a CSV file correctly formatted for the R script (Figure 2B). Specifically, the top row has the time series, and the first column contains individual sample positions on a 96-well plate. An example of such an input CSV file can be found in the supplemental material, ROS_flg22.csv or ROS_elf26.csv.

Step 3: Name the samples and/or treatments. Please note this analysis is used for data obtained from 96-well settings with time course. The experiments can be designed as 8 replicates per treatment for a total of up to 12 treatments per plate, or as 12 replicates per treatment for a total of up to 8 treatments per plate based on the rows or columns of the 96-well plate. It is important to enter the correct number of samples, either 8 or 12, because the sample number counts. If the sample count is not at 8 or 12, the code will not run. However, if there are rows of empty wells, just list them as such in the space for treatment names, e.g. naming as empty 1, empty 2... It is worth noting that for sample names, one should avoid using backslashes or similar symbols as they may cause file name issues. Additionally, one should also avoid using "NA" as a label when selecting to use ANOVA in User Input II. Unlike the LUC_2025.R script, treatments with the same name do not combine into a single group in this ROS_2025.R script. 

## User Input II: Options 1-4 can be changed per the user's specific needs.  

Option 1: Using the ANOVA test with the Tukey HSD to compare treatment luminescence sums.

Option 2: Using a two-sided t-test to compare the data. This should only be used to compare two treatments at once. The p-values shown are not adjusted to account for multiple comparisons with the same line.

Option 3: Graphing the fluorescence curves and a bar plot comparing the total sum of fluorescence.

Option 3 Addition: Adding standard deviation or standard error of mean to the graphs.

Option 4: Depending on how the plate reader records wells, the user may change the way the input is read. There are two options, one for the standard data listing of wells by A1, A2, A3... and one for wells listed by A1, B1, C1…
