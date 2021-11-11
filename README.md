# PDGMN
A Novel Feedback Vertex Set Control-based Method for Identifying Personalized Cancer Driver Genes by Using Multiplex Biomolecular Networks


## Hardware configuration
```
RAM: 16 GB
```
The code requires a good amount of RAM.


## Software requirements
- Python 3.8
- pandas 1.2.4
- numpy 1.20.1
- scikit-learn 0.24.1
- gurobi 9.1.2 
- scipy 1.6.2

### One of the most important things before running PDGMN is to make sure the "gurobi" software and its tool package "gurobipy" have been installed. (https://www.gurobi.com/)


## Usage Example

### Take BRCA cancer dataset as an example 

- Download script files of "PDGMN.py" and "utils.py", and the zip data file "data.7z" to the main directory of your project.
- Download the folder "results" and its subfolders to the main directory of your project, or you can create such folders by yourself.
- Unzip the "data.7z" to the current directory
- Set the variable cancer='BRCA' in "PDGMN.py" to perform PDGMN on BRCA cancer dataset. Note that, the variables of "netx_name" and "nety_name" indicate the names of reference networks, and these two variables are consist with the subfolders name in the "results" folder. If you want to perform PDGMN with some other reference networks, you should be aware of this.
- After that, users can directly run the command "python PDGMN.py" in the terminal to perform PDGMN for identifying personalized cancer driver genes of BRCA cancer dataset, or directly run the "PDGMN.py" file in the Pycharm IDE.
- The personalized driver genes for each patient will be sorted as a file entitled with corresponding patient ID, such as "TCGA-CH-5772-01.txt", in the folder "./results/SSMN_PPIcheng_regnetwork_BRCA/driver/".

We have provided all the necessary data files for three cancer datasets, including BRCA, LUAD and PRAD, and five different biomolecular networks as discussed in our paper. 

## Data organization

The directory of "data" contains the following directories

```
data\cancer\: contains three subfolders "data\cancer\BRCA" "data\cancer\LUAD" and "data\cancer\PRAD". Each folder sorts the corresponding data files. They are RNA-seq data file "RNA-seq.txt", somatic mutation data file "SomaticMutation.txt", and the patient IDs file "sample_lst.txt".

data\network\: contains the five different biomolecular networks.
```

The directory of "results" contains the following directories

```
results\SSMN_PPIcheng_regnetwork_BRCA: contains three subfolders "driver", "PPIcheng" and "regnetwork", in which the "driver" folder sorts the personalized driver gene files, the "PPIcheng" and "regnetwork" folders sort the sample-specific network files.

results\SSMN_PPIcheng_regnetwork_LUAD: similar to SSMN_PPIcheng_regnetwork_BRCA.
 
results\SSMN_PPIcheng_regnetwork_PRAD: similar to SSMN_PPIcheng_regnetwork_BRCA.
```
