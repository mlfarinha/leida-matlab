# Leading Eigenvector Dynamics Analysis Toolbox - Matlab

1. [Installation](#installation-of-leida-toolbox)
2. [LEiDA Start](#leida-start)
3. [LEiDA Analysis K](#leida-analysis-k)
4. [LEiDA Analysis Centroid](#leida-analysis-centroid)
5. [LEiDA Transitions K](#leida-transitions-k)
6. [LEiDA State Time Courses](#leida-state-time-courses)

---

## Installation of LEiDA Toolbox

Start by creating a new folder where all the LEiDA code and analyses will be saved. Then to install the LEiDA Toolbox unzip the folder LEiDA_Functions to the created folder and add that folder (with all its functions) to your Matlab path. In your Matlab command terminal type the following:

`addpath(genpath('my_LEiDA_directory\LEiDA_Functions'))`

If you would like to save the LEiDA_Functions folder to your Matlab path for later use:

1. click *HOME* in your Matlab dashboard;
2. click *Set Path*;
3. click *Save* (making sure you see the directory my_LEiDA_directory\LEiDA_Functions in your Matlab search path pop-up window);
4. click *Close*;

> Note that the LEiDA Toolbox requires the Statistics and Machine Learning Matlab Toolbox.

---

## LEiDA Start

The function ***LEiDA_Start.m*** runs the first stage of the LEiDA analyses. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***Data_directory***: directory of the folder with the parcellated fMRI data; this folder does not need to be inside the folder *LEiDA_directory*;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;
- ***Conditions_tag***: tag of the conditions given in the imaging files; it is expected that the file names of the parcellated data files contain a tag stating the condition to which a given participant belongs to; the tags should correspond exactly to the tags used in the file names; there can be any number of conditions;
- ***Parcellation***: parcellation template applied to the segment the fMRI data. Available options include 'AAL116', 'AAL120', 'dbs80' and 'glasser360';
- ***N_areas***: number of brain areas of the parcellation to consider for analysis;
- ***TR***: TR of the fMRI data;
- ***Tmax***: maximum number of TRs considering all fMRI sessions;
- ***apply_filter***: apply temporal filtering to the parcellated data (0: no; 1: yes);
- ***flp***: lowpass frequency of filter (default 0.1);
- ***fhi***: highpass frequency of filter (default 0.01);
- ***Paired_tests***: experimental paradigm (0: subjects in different conditions are not the same; 1: subjects are the same across conditions);
- ***CortexDirection***: direction to plot the PL states/brain ('SideView' or 'TopView').

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Creates a directory to store the results from the LEiDA analyses;
- Computes the leading eigenvectors of the data and saves them in the file *LEiDA_EigenVectors.mat*;
- Clusters the leading eigenvectors of all participants and saves the clusters in the file *LEiDA_Clusters.mat*;
- Computes the fractional occupancy, performs hypothesis (permutation) tests to compare the mean fractional occupancy across conditions and saves the results in the file *LEiDA_Stats_FracOccup.mat*;
- Plots figure with two-sided p-values obtained from the comparison of the mean fractional occupancy, plots figures with barplots of the fractional occupancy values and saves the figures in png and fig format;
- Computes the dwell time, performs hypothesis (permutation) tests to compare the mean dwell time across conditions and saves the results in the file *LEiDA_Stats_DwellTime.mat*;
- Plots figure with two-sided p-values obtained from the comparison of the mean dwell time, plots figures with barplots of the dwell time values and saves the figures in png and fig format;
- Plots the PL states (centroids) obtained by clustering the leading eigenvectors and saves the figure in png and fig format;

---

## Format of the inputs

The input data should be contained in the folder specified in the variable *Data_directory*. This folder should contain the time series data for each participant as a separate file. The file name of each participant should contain a tag specifying the condition of each participant to allow the computation of intergroup differences. The time series data can be supplied in the following formats:

- A (no. of brain areas $\times$ no. of time points) matrix with double or single precision numeric elements for each subject. The file can be either a Matlab file (.mat), a text file (.txt) or a .1D file. If the files are Matlab structure arrays, the time series fMRI data must be enclosed in a field named data.

---

## LEiDA Analysis K

The function ***LEiDA_AnalysisK.m*** plots a number of relevant figures to analyse the LEiDA results for a specific number of clusters. Based on the based on the outputs from the previous analyses (*LEiDA_Start.m*) the user should select one or more values for K to explore with greater detail. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states to consider for further analyses;
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;
- ***Parcellation***: parcellation template applied to the segment the fMRI data;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Plots several plots of the selected K:
    - PL states as links in cortex between the brain areas with positive values in the centroid vector;
    - PL states as nodes in cortex, where each brain area is coloured and scaled according to its value in the centroid vector;
    - PL states rendered in a 3D glass brain indicating the overlap with the functional brain networks defined by Yeo et al., (2011);
    - PL states in matrix format;
    - PL states in vector format with numbered areas;  
    - Boxplot of the fractional occupancy values of all K PL states;
    - Boxplot of the dwell time values of all K PL states;
    - PL states rendered in 3D glass brain, in vector format and boxplots of the fractional occupancy and dwell time values;
    - Overlap between the K PL states and the 7 RSNs defined by Yeo et al., (2011);

The plots are saved in new directory created according to the value defined for K.

---

## LEiDA Analysis Centroid

The function ***LEiDA_AnalysisCentroid.m*** should be used to analyse a PL state for a given value of K selected by the user. This function plots a number of relevant figures to analyse this PL state. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states (K) to consider for further analyses;
- ***Centroid***: define the PL state to consider for further analyses (1 <= Centroid <= SelectK);
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;
- ***Parcellation***: parcellation template applied to the segment the fMRI data;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Plots several plots of the selected PL state:
    - PL state in vector format with label of each brain area according to the specified parcellation;
    - PL state in vector format with brain areas organised by contribution;
    - Boxplot of the fractional occupancy values for the selected PL state;
    - Boxplot of the dwell time values for the selected PL state;
    - PL state in vector format, rendered in a 3D glass brain, in matrix format and boxplots of the fractional occupancy and dwell time values;

The plots are saved in new directory created according to the value of K and the specified PL state.

---

## LEiDA Transitions K

The function ***LEiDA_TransitionsK.m*** should be used to analyse the state-to-state transitions considering a state space of dimension to be defined by the user, i.e., by the value selected for K. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states (K) to consider for further analyses;
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Computes the Transition Probability Matrix (TPM) for each participant and performs permutation tests to compre the state-to-state transition probabilities between conditions. The results are saved in the file *LEiDA_Stats_TransitionMatrix.mat* to the folder corresponding to selected K;
- Plots the mean TPM for each condition;
- Plots the intergroup differences between state-to-state transition probabilties;

---

## LEiDA State Time Courses

The function ***LEiDA_TransitionsK.m*** should be used to analyse the state time courses defined for a specified value of K across participants and for a given subject. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states (K) to consider for further analyses;
- ***Subject***: define the subject whose state time course should be plotted. The user should provide the filename or an ID/number that uniquely identified the subject;
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Plots the state time courses considering K PL states for all participants and separates them by condition;
- Plots the state time course for a given subject as cluster blocks;
- Plots the state time course for a given subject as stairs plot;