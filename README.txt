ClonePhyTester_v0.1.1
Updated August 16, 2019
==================
ClonePhyTester was developed by Sudhir Kumar [1].

ClonePhyTester computes error rates (MLTED [2], RF [3], TreeVec [4], and mutation order error rates; see ref. [1] for the detail of these scores) for given clone phylogenies, inferred from multi-sample bulk sequencing data obtained from a single patient. These error rates can be compared with the clone prediction methods evaluated in ref. [1] (CloneFinder, MACHINA, Treeomics, LICHeE, MixPhylogeny, PhyloWGS, and Cloe). 

ClonePhyTester is built for Windows 10. Linux version is currently not available. 



Dependencies
==================
1. R (version 3.6.0, 3.5.3, and 3.5.2 were tested)
    R dependencies:
    ape, maps, phytools, ade4, phangorn, jsonlite, googleVis, treespace, ggplot2, reshape2, and gdata

Note: To install these R packages, you may run the following in R:   install.packages(c("ape","maps","phytools","ade4","phangorn","jsonlite","googleVis","treespace","ggplot2","reshape2","gdata"))

    If gdata is not properly installed and the following error message is encountered, “Unable to load perl libraries needed by read.xls() to support XLX files, ” you may need to download Perl (version 5.12.30) from https://github.com/dwimperl/perl-5.12.3.0. Then, unzip into Windows(C:) folder. You can test if all dependencies were correctly installed by running the following command in R:

library(ape) 
library(maps)
library(phytools)
library(ade4)
library(phangorn)
library(jsonlite)
library(googleVis)
library(treespace)
library(ggplot2)
library(reshape2)
library(gdata)

2. python (version 2.7.13)
    python  dependencies:
    Biopython, Pandas, numpy, and ete2

Note: If the installation of these python packages is not easy, you may want to use Anaconda (4.4.0, 64-bit) for python 2.7 (https://www.anaconda.com/distribution/). If you use Anaconda, run the following commands from command prompt.

conda install -c anaconda biopython
conda install -c anaconda pandas
conda install -c anaconda numpy
conda install -c etetoolkit ete2 

Or, the following commands can also work.

pip install Biopython
pip install pandas
pip install numpy
pip install ete2

3. MLTED (https://github.com/khaled-rahman/MLTED)[2]
By using a command prompt, run the following command in the ClonePhyTester directory:
g++ -std=c++11 main.cpp -o main

Note: main.exe should be created in the ClonePhyTester directory. Prior installation of g++ is required. You can use MinGW (http://www.mingw.org/). For the details of the installation steps, this website may be helpful (http://www.codebind.com/cprogramming/install-mingw-windows-10-gcc/).

4. MEGA-CC 
Please download the latest version from https://www.megasoftware.net/.



How to use 
==================
1. Select a clone prediction method (which is not included in ref. [1]) and infer clone sequences on given read count information.
You need to select a clone prediction method that is designed for the analysis of bulk sequencing data of multiple samples from a single patient (e.g., PhyloWGS). The input information should be read count, observed SNV frequencies, or cancer cell fraction. Simulated read count data are given in the directory of G7, G12, MA, and P10. These data was used in ref. [1] and the details are described in the section of ‘Simulated datasets.’

2. Save your inferred clone sequences in the MEGA format (https://www.megasoftware.net/). 
Note that 'A' and 'T' denote wild-type and mutant-type base assignments, respectively. The order of SNV should be the same as the read count data. Do not include normal (germline) sequences.

For example: 
Let us assume that you inferred some clone sequences by using a few read count data in G7, G12, P10, and MA directory. We saved these inferred clone sequences in the MEGA format in the directory of “Example.” 

3. List the mega files that contain your inferred clone sequences.
Open “sequenceList.txt” in the ClonePhyTester directory. Do not change the first two columns (i.e., “Sim” and “Data”). The “Sim” section indicates which simulation was used to generate a dataset, i.e., G7, G12, P10, or MA. Under the “Data” column, the ID of dataset was listed. List your inferred clone sequences in third column (“Inferred Seq”) accordingly. 

For example: 
The example inferred clone sequences were listed in Test_sequenceList.txt, which can be found in the Example directory. 

4. Make a control file for ClonePhyTester. 
List the paths to the file that lists your clone sequences, your Rscript.exe, and main.exe (MLTED).
For example: 
For the example datasets, we created test.input that can be found in the Example directory. If you like to test these example datasets, then you should change the paths accordingly. 

5. Run clonephytester.py.
Open anaconda prompt or command prompt and run the command, 
python clonephytester.py [your input file]

For example: 
To perform the example data analysis, try:
python clonephytester.py Example/test.input

Note: ClonePhyTester may not run from Dropbox or Googledrive.



Output file
==================
The output files are found in the folder that contains the control file (e.g., test.input). 
1. Summary plot (jpeg file). 
All scores (MLTED, Mutation Orders, TreeVac, and RF) are plotted with the other methods that were tested in ref. [1], i.e., CloneFinder, MACHINA, TreeOmics, LICHeE, MixPhy, PhyloWGS). The scores of your inferred clone sequences will be shown under “YourMethod” in the plot. If you do not use all the input files that are found in G7, G12, P10, and MA directories, some methods and some panels may not show plots in the figure. We use only datasets for which inferred clone sequences are provided. Also, when inferred clone sequences do not produce a clone phylogeny (e.g., star-like phylogeny), that dataset will be excluded. If your inferred clone sequences show higher error rates than the best performing methods in ref [1] (i.e., CloneFinder, MACHINA, Treeomics, and LICHeE), we recommend to not use that method. 

2. Summary table (excel file). 
The error rates for each dataset are provided in an excel file, together with the other methods. 
For the example data analysis, it will produce test_summary.jpeg and test_summary.xlsx.



Simulated datasets
==================
Bulk sequencing data was simulated (read count).
Four different evolutionary scenarios were used (G7, G12, P10, and MA).
Correct clone sequences are given in the MEGA format and can be open by using the MEGA software (https://www.megasoftware.net/).
'A' and 'T' denote wile-type and mutant-type base assignments, respectively.
The order of SNV is the same between the bulk sequencing data and correct clone sequences

Graphical visualization of the results in ref. [1]
==================================================
If you like to visualize the results of clone prediction methods tested in ref. [1] without your method, please use Plots-withoutYourMethod.r. To run this code, you need to change lines 5 and 6, which assign the path to ClonePhyTester directory and output file name.

Reference:
[1] Miura S, Vu T, Deng J, Buturla T, Choi J, Kumar S: Power and pitfalls of computational methods for inferring clone phylogenies and mutation orders from bulk sequencing data. bioRxiv 2019:697318. 
[2] Karpov N, Malikic S, Rahman K, Sahinalp S.C.: A Multi-labeled Tree Edit Distance for Comparing "Clonal Trees" of Tumor Progression. In: 18th International Workshop on Algorithms in Bioinformatics (WABI 2018). Volume 113: pp 22:1--22:19 
[3] Wen D, Yu Y, Zhu J, Nakhleh L: Inferring Phylogenetic Networks Using PhyloNet. Syst 747 Biol 2018, 67(4):735-740. 
[4] Kendall M, Eldholm V, Colijn C: Comparing phylogenetic trees according to tip label categories. bioRxiv 2018:251710.
