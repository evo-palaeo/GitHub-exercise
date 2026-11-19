### **Overview**

This repository contains Python code for analysing biological sequences using Biopython. It is intended to be used in a practical to demonstrate core Git and GitHub skills, including version control, branching, and working with remote repositories. Students modify and extend the scripts to fix bugs, add new features, and practise good coding workflows.

### **Repository Structure**

| File          | Description                                           |
|---------------|--------------------------------------------------------|
| `translate.py` | Main script for ORF detection and translation (contains intentional bugs ⚠️) |
| `input.fas`    | Small synthetic FASTA file for trying the script      |
| `README.md`    | Project documentation                                |


### **Features**

- Reads and processes DNA sequences from FASTA files
- Identifies open reading frames
- Translates sequences into proteins
- Computes sequence statistics and stores them in a pandas DataFrame

### **Installation**

`pip install biopython`

`pip install pandas`

### **Student Exercise**

- **Fork a repository** <br>Navigate to https://github.com/evo-palaeo/GitHub-exercise and create a fork under your GitHub account. 
- **Clone your fork** <br>Use VS Code to clone your fork to your local computer. 
- **Fix bugs in translate.py** <br>Open translate.py in VS Code, identify and fix the bugs, then stage, commit, and push your changes to your fork. 
- **Create a development branch** <br>Create a new branch with a descriptive name for your development work. 
- **Modify the script in the development branch** <br>Add a Python `input()` function to ask the user for a FASTA filename. After making the change, stage, commit, and push it to the development branch. 
- **Update the main branch** <br>In the main branch, add Python code to plot Sequence Length vs Molecular Weight using matplotlib. 
- **Merge branches** <br>Merge the development branch into the main branch, resolving any conflicts if necessary.
