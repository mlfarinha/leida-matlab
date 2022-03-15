# Leading Eigenvector Dynamics Analysis Toolbox - Matlab

1. [Installation](#installation-of-leida-toolbox)
2. [LEiDA Master](#leida-master)
3. [LEiDA Optimal K](#leida-optimal-k)

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

## LEiDA Master

The function ***LEiDA_Master.m*** runs the first stage of the LEiDA analyses. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***Data_directory***: directory of the folder with the parcellated fMRI data; this folder does not need to be inside the folder *LEiDA_directory*;
- ***Conditions_tag***: tag of the conditions given in the imaging files; it is expected that the file names of the parcellated data files contain a tag stating the condition to which a given participant belongs to; the tags should correspond exactly to the tags used in the file names; there can be any number of conditions;
- ***Parcellation***: parcellation template applied to the segment the fMRI data;

> This version is only working with the *AAL116* parcellation at the moment; this will be updated soon;

- ***N_areas***: number of brain areas of the parcellation to consider for analysis;
- ***TR***: TR of the fMRI data;
- ***Tmax***: maximum number of TRs considering all fMRI sessions;
- ***apply_filter***: apply temporal filtering to the parcellated data (0: no; 1: yes);
- ***flp***: lowpass frequency of filter (default 0.1);
- ***fhi***: highpass frequency of filter (default 0.01);
- ***Paired_tests***: experimental paradigm (0: subjects in different conditions are not the same; 1: subjects are the same across conditions);
- ***CortexDirection***: direction to plot the FC states/brain ('SideView' or 'TopView').

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Creates a directory to store the results from LEiDA analyses;
- Computes the leading eigenvectors of the data and saves them in the file *LEiDA_EigenVectors.mat*;
- Clusters the leading eigenvectors of all participants and saves the clusters in the file *LEiDA_Clusters.mat*;
- Computes the fractional occupancy, performs hypothesis tests to compare the mean fractional occupancy across conditions and saves the results in the file *LEiDA_Stats_FracOccup.mat*;
- Plots figure with two-sided p-values obtained from the comparison of the mean fractional occupancy, plots figures with barplots of the fractional occupancy values and saves the figures in png and fig format;
- Computes the dwell time, performs hypothesis tests to compare the mean dwell time across conditions and saves the results in the file *LEiDA_Stats_DwellTime.mat*;
- Plots figure with two-sided p-values obtained from the comparison of the mean dwell time, plots figures with barplots of the dwell time values and saves the figures in png and fig format;
- Plots the centroids obtained by clustering the leading eigenvectors, computes their overlap with the networks defined by Yeo et al., (2011) and saves the figure in png and fig format;

---

## LEiDA Optimal K

The function ***LEiDA_OptimalK.m*** runs the analysis for the optimal number of FC states selected based on the outputs from the previous analyses (*LEiDA_Master.m*). This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***OptimalK***: define the number of FC states to consider for further analyses;
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***Parcellation***: parcellation template applied to the segment the fMRI data; this version is only working with the *AAL116* parcellation at the moment; this will be updated soon;
- ***N_areas***: number of brain areas of the parcellation to consider for analysis;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Creates a directory to store the results for the selected optimal K;
- Plots several summary plots of the optimal K:
    - FC states rendered in a 3D glass brain indicating the overlap with the Yeo functional brain networks;
    - FC states in matrix format;
    - FC states in vector format with labels of parcellation;
    
    > At this stage, only the labels of 90 areas of the AAL parcellation are provided. This will be updated to include other parcellations and the corresponding labels.
    
    - FC states in vector format with numbered parcels;
    - FC states in vector format highlighting the contribution of each brain area;
    - Boxplot of the fractional occupancy values of all FC states;
    - Boxplot of the fractional occupancy values of each FC state;
    - Boxplot of the dwell time values of all FC states;
    - Boxplot of the dwell time values of each FC state;
    - Vector, 3D glass brain, matrix and boxplots of the fractional occupancy and dwell time values for each FC state;
    - 3D glass brain, vector and boxplots of the fractional occupancy and dwell time values for all FC states;
- Computes the transition probability matrix, performs hypothesis tests to compare the mean state-to-state transition probabilities across conditions and saves the results in the file *LEiDA_Stats_TransitionMatrix.mat* to the folder corresponding to selected K;
- Plots the mean transition probability matrix and the differences in state-to-state transition probabilities between conditions.