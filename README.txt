ClonePhyTester_v0.1.1  
(Copyright 2019, Authors and Temple University; see license below)

Updated October 28, 2019
==================

The ClonePhyTester pipeline has been developed by Sayaka Miura. It is written in Python3.7 and Snakemake and work smoothly on Linux/Unix at present. You are free to download, modify, and expand this pipeline under a permissive license similar to the BSD 2-Clause License (see below). ClonePhyTester has a few standard dependencies that are noted below.

The ClonePhyTester pipeline was developed for conducting a benchmarking analysis (Miura et al., Genome Biology, in review; ref. 1) in which we evaluated the performance of seven computational methods (CloneFinder, MACHINA, Treeomics, LICHeE, MixPhylogeny, PhyloWGS, and Cloe) in producing clone phylogenies for many simulated datasets. We assessed the accuracy of tested methods metrics in determining the order of mutations and the branching pattern within the reconstructed clone phylogenies (four different performance metrics are used: MLTED, RF, TreeVec, and mutation order error rates). For scientific rationale, details, and references, please consult Miura et al.(Genome Biology, in review; ref. 1). 

In brief, if you download ClonePhyTester, you will automatically receive all the simulated datasets and primary results produced by the seven methods we tested. That is, you will not need to download and collect the datasets from many different original sources and programs.  If you execute the ClonePhyTester pipeline as is after downloading, it will replicate all the results, because the results from all the seven methods applied to all the simulated datasets are also packaged in the download.  ClonePhyTester contains all the utilities to compute the performance metrics, and it will produce graphical visualizations and summaries. For all other uses, see below for How to use ClonePhyTester.

Dependencies
==================
1. R (version 3.6.1 was tested)
    R dependencies:
    Rcpp, ape, maps, curl, phytools, ade4, phangorn, treespace, ggplot2, reshape2, and gdata
Note: To install these R packages, you may run the following in R:   install.packages(c("Rcpp","ape","maps","curl","phytools","ade4","phangorn","treespace","ggplot2","reshape2","gdata"))

While downloading the packages, if "curl" or "phytools" do not install properly, the dependencies, libcurl and magick, might need to be installed using the following commands in the terminal line before the packages get installed again:
	sudo apt-get install libcurl4-openssl-dev
	sudo apt-get install libmagick++-dev
	
While downloading the packages, if "treespace" do not install properly, the dependencies, unit, libgdal, and libudunits2, might need to be installed using the following commands in the terminal line before the packages get installed again:
	sudo apt-get install libssl-dev
	sudo apt-get install libudunits2-dev
	sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
	sudo apt-get install libgdal-dev

You can test if all dependencies were correctly installed by running the following commands in R:
library(Rcpp)
library(ape) 
library(maps)
library(curl)
library(phytools)
library(ade4)
library(phangorn)
library(treespace)
library(ggplot2)
library(reshape2)
library(gdata)

2. python (version 3.7)
    python  dependencies:
    Snakemake, Biopython, Pandas, numpy, xlrd, xlsxwriter, and ete3
Note: If the installation of these python packages is not easy, you may want to use Anaconda for python 3.7 (https://www.anaconda.com/distribution/). Or, you can try python3-pip, by using the following command:
sudo apt install python3-pip
sudo python3 -m pip install snakemake
sudo python3 -m pip install Biopython
sudo python3 -m pip install pandas
sudo python3 -m pip install numpy
sudo python3 -m pip install ete3
sudo python3 -m pip install xlrd
sudo python3 -m pip install xlsxwriter

3. MLTED (https://github.com/khaled-rahman/MLTED)[2]
By using a command prompt, run the following command in the ClonePhyTester directory:
g++ -std=c++11 main.cpp -o main

4. MEGA-CC 
Please download MEGA-CC from,
Windows: https://www.megasoftware.net/releases/MEGA_X_CC_10.1.7_win64_setup.exe
Ubuntu/Debian: https://www.megasoftware.net/releases/megax-cc_10.1.7-1_amd64.deb 
Redhat/Fedora: https://www.megasoftware.net/releases/megax-cc-10.1.7-1.x86_64.rpm
macOS: https://www.megasoftware.net/releases/MEGAX_CC_10.1.7_installer.pkg

5. Java
Please download the latest version of Java using the command prompt.

If you wish to use ClonePhyTester to compare the performance of another method (new or existing) with the currently analyzed seven methods for the simulated dataset we employed, then ClonePhyTester will need you to give the clone sequences that are produced by the new method. With that information, ClonePhyTester will calculate all the performance metrics and produce graphical visualizations and summaries that will compare the new method with the existing seven methods. For adding new metrics or datasets, you should be able to easily modify ClonePhyTester, as it is ClonePhyTester is written in python.

How to use 
==================
1. Select a clone prediction method (which is not included in ref. [1]) and infer clone sequences on given read count information.
You need to select a clone prediction method that is designed for the analysis of bulk sequencing data of multiple samples from a single patient (e.g., PhyloWGS). The input information should be read count, observed SNV frequencies, or cancer cell fraction. Simulated read count data are given in the directory of G7, G12, MA, P10, G7cna, MA50, TGlinear, TGstep, and TGconst. These data were used in ref. [1] and the details are described in the section of "Simulated datasets."

2. Save your inferred clone sequences in the MEGA format (https://www.megasoftware.net/). 
Note that 'A' and 'T' denote wild-type and mutant-type base assignments, respectively. The order of SNV should be the same as the read count data. Do not include normal (germline) sequences.
For example: 
Let us assume that you inferred some clone sequences by using a few read count data in G7, G12, P10, and MA directory. We saved these inferred clone sequences in the MEGA format in the directory of "Example/Sequence_List." 

3. List the mega files that contain your inferred clone sequences.
Open "SequenceList.txt" in the ClonePhyTester directory. Do not change the first two columns (i.e., "Sim" and "Data"). The "Sim" section indicates which simulation was used to generate a dataset, e.g., G7 and MA. Under the "Data" column, the ID of dataset was listed. List your inferred clone sequences in third column ("Inferred Seq") accordingly. 

For example: 
The example inferred clone sequences were listed in Test_sequenceList.txt, which can be found in the Example directory. 

4. Make a control file for ClonePhyTester. 
List the paths to the file that lists your clone sequences, your Rscript, and ./main (MLTED).
For example: 
For the example datasets, we created Test.input that can be found in the Example directory. If you like to test these example datasets, then you should change the paths accordingly. 

5. Run clonephytester.py.
Open terminal and run the command, 
python clonephytester.py [your input file]

For example: 
To perform the example data analysis, try:
python clonephytester.py Example/Test.input

If you like to run ClonePhyTester with Snakemake, Open terminal and run the command, 
snakemake --config filename=[Your input file] 

For example: 
To perform the example data analysis, try:
snakemake --config filename=Example/test.input

To run in parallel, add -j [thread_count] to the end of the command.

Note: ClonePhyTester may not run from Dropbox or Googledrive.

Output file
==================
The output files are found in the folder that contains the control file (e.g., test.input). 
1. Summary plot (jpeg file). 
All scores (MLTED, Mutation Orders, TreeVac, and RF) are plotted with the other methods that were tested in ref. [1], i.e., CloneFinder, MACHINA, TreeOmics, LICHeE, MixPhy, PhyloWGS). The scores of your inferred clone sequences will be shown under "YourMethod" in the plot. If you do not use all the input files that are found in G7, G12, P10, MA, MA50, G7cna, TGlinear, TGstep, and TGconst directories, some methods and some panels may not show plots in the figure. We use only datasets for which inferred clone sequences are provided. Also, when inferred clone sequences do not produce a clone phylogeny (e.g., star-like phylogeny), that dataset will be excluded. If your inferred clone sequences show higher error rates than the best performing methods in ref [1] (i.e., CloneFinder, MACHINA, Treeomics, and LICHeE), we recommend to not use that method. Note that G7cna contains CNAs, which may lose SNVs, so MLTED and error rates of ordering mutations were not computed. Also, methods tested in ref [1] substantially underestimated clone count for TGlinear, TGstep, and TGconst; we did not compete RF and TreeVec. 

2. Summary table (excel file). 
The error rates for each dataset are provided in an excel file, together with the other methods. 
For the example data analysis, it will produce test_summary.jpeg and test_summary.xlsx.

Simulated datasets
==================
Bulk sequencing data was simulated (read count).
Seven different evolutionary scenarios were used (G7, G12, P10, MA, TG, G7cna, and MA50). G7cna datasets are G7 datasets with CNAs, and MA50 datasets are MA datasets with the sequencing read depth of 50.
Correct clone sequences are given in the MEGA format and can be open by using the MEGA software (https://www.megasoftware.net/).
'A' and 'T' denote wile-type and mutant-type base assignments, respectively.
The order of SNV is the same between the bulk sequencing data and correct clone sequences

Tips for reuse of clonephytester.py
===================================
The clonephytester.py program is designed to compute scores described in ref. [1] by using the same datasets as in ref [1]. Since clonephytester.py is written in Python3, you can easily modify the code for your purpose. Below, we provide two examples. 

Example A: Using datasets that are generated by users
By making a few changes in clonephytester.py, the same scores as in ref [1] can be computed for datasets generated by any user. See below for instructions for a set of example datasets and input files that are included in the Example/Tips_for_modification_of_clonephytester/Example_A-user_generated_data directory. Inferred clone genotypes (sequences) are InfA1.meg, InfA2.meg, and InfA3.meg, which can be found in SimA/Inferred directory. The format of the inferred clone sequences is the same as the original clonephytester.py, i.e., "A" and "T" denote wild-type base and mutant bases, respectively, and the order of SNVs need be the same in all the inferred clone sequences. 

Step 1: Prepare files with true clone genotypes
True clone genotypes need to be saved in a MEGA format. In the alignment, the order of SNVs needs to be the same as that in the inferred clone sequences, and wild-type and mutant bases need to be indicated by "A" and "T", respectively. For our example, we deposited the true clone genotypes in the SimA/True directory, i.e., A1_True.meg, A2_True.meg, and A3_True.meg that correspond to their inferred clones, InfA1.meg, InfA2.meg, and InfA3.meg, respectively found in the SimA/Inferred directory.  

Step 2: Prepare the sequence list file
The sequence list file contains information on simulation ID, the location of true sequences (A1_True.meg, A2_True.meg, and A3_True.meg for this example), and the location of inferred sequences (InfA1.meg, InfA2.meg, and InfA3.meg for this example). Any name can be assigned for the simulation ID, but it should be short. A sequence list file should look as follows,
Sim	Data	Inferred Seq
A	True1.meg	Inf1.meg
A	True2.meg	Inf2.meg
B	TrueB.meg	InfB.meg 
For our example, the sequence list file is SimA_SequenceList.txt and found in the directory (Example/Tips_for_modification_of_clonephytester/Example_A-user_generated_data).

Step 3: Prepare a control file (.input)
The format of a control file is the same as for the original clonphytester.py. See SimA.input for our example.

Step 4: Prepare a template of summary excel file (output excel file)
A template excel file needs to contain six sheets that are named as "MLTED", "RF", "TreeVec, "Sequential_MutPair_count", "parallel_MutPair_count", and "concurrent_MutPair_count", respectively. Alternatively, a user can modify Summary.xlsx in the main directory of ClonePhyTester. Each sheet contains the information of Simulation ID and Data ID (true clone sequence file), which need to be identical to those listed in the sequence list file (SimA_SequenceList.txt for this example). See Summary_example_A.xlsx for this example. 

Step 5: Edit clonephytester.py code
(1) Assign your simulation ID (line 21 of clonephytester.py)
The dictionary, "Sim2TrueFile", connects each true clone sequence file to its inferred clone sequence file by reading the sequence list file. The key-item pair for this dictionary is the simulation ID and the word "ID". The simulation ID should be the same as that in the sequence list file. For our example, the new command in the modified clonephytester.py program is,
Sim2TrueFile={"SimA":"ID"}
 
(2) Change the template summary excel file name (Line 104 in clonephytester.py)
For our example, the new command is,
xl = pd.ExcelFile('Example/Tips_for_modification_of_clonephytester/Example_A-user_generated_data/Summary_example_A.xlsx')
These edits have been included in ExampleA-clonephytester-user_generated_data.py for your reference. 

Example B: Calculating scores defined by users
The ClonePhyTester framework can be easily modified if one wishes to add or delete score computations. In this case, the preparation of true clones, inferred clones, sequence list, and control files are the same as in the original clonephytester.py. If a user would like to use their own datasets, the same procedure that is described in Example A is necessary. We explain the steps for defining new scores. In this example, we assume that a new score is clone count error, i.e., the difference of inferred clone count from its true clone count. We call this score "Clonecount."
 
Step 1: Prepare a template of summary excel file (output excel file)
The content of each sheet is the same as the other template of summary excel file (see Summary.xlsx as an example). For a score defined by a user, a new sheet needs to be created in the same excel file, and the name of the sheet is the name of score. For the example of "Clonecount" score, we create a sheet named "Clonecount." In this example, we assume we do not need the other scores, so our excel file contains only this sheet. See Summary_example_B.xlsx as a reference (at Example/Tips_for_modification_of_clonephytester/Example_B-user_specified_score). 

Step 2: Edit clonephytester.py code
(1) Add a code to compute scores defined by users
In the original clonephytester.py, Lines 69 to 96 of this script are the commands to compute scores used in ref [1], i.e., RF, MLTED, TreeVec, and mutation error rate. If users do not need these scores, just delete the lines containing those names. In this example, we decided not to compute other scores so we deleted all of these lines from the original clonephytester.py file. 
Here, all the scores are stored in the dictionary, "ResAll," in which the key is the name of score ("Clonecount" for this example) and items are the score for an inferred clones. To calculate a new score, users can use "functions" in Python and can be saved in Functions.py. In this example, we added the following code in the Functions.py:
###Example of adding your own score computation###
def Compute_clone_count_error(TruMegFile,InfMegFile):
    TcloLs, Tcloseq = ReadMegSeq(TruMegFile)
    IcloLs, Icloseq = ReadMegSeq(InfMegFile)
    CCerror=len(IcloLs)-len(TcloLs)
    return CCerror
##################################################

And the following line was added into the original clonephytester.py at line 92-94:
####Count the number of inferred clones and calculate the difference from the true count####
CloneCount_error=Functions.Compute_clone_count_error(TruMegTMP,InfMegTMP)
ResAll['CloneCount']=CloneCount_error

(2) Change the template summary excel file name (Line 107 in clonephytester.py)
For our example, the new line is,
xl = pd.ExcelFile('Example/Tips_for_modification_of_clonephytester/Example_A-user_generated_data/Summary_example_B.xlsx')
These edits were saved in ExampleB-clonephytester-user_specified_score.py for a reference. 


Reference:
[1] Miura S, Vu T, Deng J, Buturla T, Choi J, Kumar S: Power and pitfalls of computational methods for inferring clone phylogenies and mutation orders from bulk sequencing data. bioRxiv 2019:697318. 
[2] Karpov N, Malikic S, Rahman K, Sahinalp S.C.: A Multi-labeled Tree Edit Distance for Comparing "Clonal Trees" of Tumor Progression. In: 18th International Workshop on Algorithms in Bioinformatics (WABI 2018). Volume 113: pp 22:1--22:19 
[3] Wen D, Yu Y, Zhu J, Nakhleh L: Inferring Phylogenetic Networks Using PhyloNet. Syst 747 Biol 2018, 67(4):735-740. 
[4] Kendall M, Eldholm V, Colijn C: Comparing phylogenetic trees according to tip label categories. bioRxiv 2018:251710.

--------
Copyright 2019, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
