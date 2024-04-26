# GS1_GenomeAnalysingProgram

1. Install Required Python Packages
Make sure to install the required Python packages using pip if you haven't already. 
pip install panel biopython bokeh numpy scikit-learn

2. Install Muscle Sequence Alignment Tool
Download and install the Muscle sequence alignment tool on your system.
You've referenced muscle3.8.31_i86win32.exe, which is the Windows executable for Muscle. Follow these steps to set up Muscle:
Download Muscle: Visit the Muscle website and download the appropriate version for your operating system.
Install Muscle: After downloading, extract the Muscle executable (muscle3.8.31_i86win32.exe for Windows) to a directory included in your system's PATH environment variable, or place it in the same directory as your Python script.

3. Prepare FASTA file for Large Sequences:
Ensure that the FASTA file contains the sequences you want to analyze.
The program is capable of handling large sequence files, although processing time may be significant depending on the size of the input.
Upload FASTA Files:
Use the program's file upload functionality to select and upload the desired FASTA file(s) containing the sequences.
Mutation Analysis:
Specifically for mutation analysis, the program requires exactly two input files.
Select and upload the two specific FASTA files that you want to use for mutation analysis.


** We use file name "genomeversion2.py" as lasted file**
If you want to see our journey please look at "genomev1.py"
