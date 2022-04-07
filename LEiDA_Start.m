function LEiDA_Start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function runs the first stage of the LEiDA analyses. It should be
% used to select the number of FC states appropriate for a more detailed
% analysis. This function contains two sections: (A) user input parameters;
% and (B) code to run LEiDA. The user should only need to adapt section A.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Run Leading Eigenvector Dynamics Analysis:
%    - Compute the leading eigenvectors for all participants
%    - Cluster the leading eigenvectors of all participants
%    - Analysis of Fractional Occupancy values
%    - Analysis of Dwell Time values
%    - Pyramid of FC states
%
% Tutorial: README.md
% Version:  V1.0, April 2022
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: USER INPUT PARAMETERS

% Directory of the LEiDA toolbox folder:
LEiDA_directory = 'D:/LEiDA_Toolbox/';
% Directory of the folder with the parcellated neuroimaging data:
Data_directory = 'D:/LEiDA_Toolbox/ABIDE_dparsf_all_aal116/';
% Name of the run to be used to create the folder to save the data:
run_name = 'ABIDE_dparsf_all_AAL116';
% Tag of conditions given in the parcellated image files:
Conditions_tag = {'Control','Autism'};
% Parcellation applied to the imaging data (see tutorial):
Parcellation = 'AAL116';
% Number of brain areas to consider for analysis;
N_areas = 90;
% Repetition time (TR) of the fMRI data (if unknown set to 1):
TR = 1;
% Maximum number of TRs for all fMRI sessions:
Tmax = 315;
% Apply temporal filtering to data (0: no; 1: yes)
apply_filter = 0;
% Lowpass frequency of filter (default 0.1):
flp = 0.1;
% Highpass frequency of filter (default 0.01):
fhi = 0.01;
% Experimental paradigm (0: subjects in different conditions are not the
% same; 1: subjects are the same across conditions)
Paired_tests = 0;
% Direction to plot the FC states/brain ('SideView' or 'TopView'):
CortexDirection = 'SideView';

% AFTER FILLING IN THE INPUT PARAMETERS AND ADDING LEiDA_Function FOLDER TO
% YOUR MATLAB PATH:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

%% B: RUN LEIDING EIGENVECTOR DYNAMICS ANALYSIS

% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)

% Create a directory to store the results from the current LEiDA run
if ~exist([LEiDA_directory 'LEiDA_Results_' run_name '/'], 'dir')
    mkdir([LEiDA_directory 'LEiDA_Results_' run_name '/']);
end
leida_res = [LEiDA_directory 'LEiDA_Results_' run_name '/'];

% Compute the leading eigenvectors of the data
LEiDA_data(Data_directory,leida_res,N_areas,Tmax,apply_filter,flp,fhi,TR);

% Cluster the leading eigenvectors of all subjects
LEiDA_cluster(leida_res);

% Compute the fractional occupancy and perform hypothesis tests
LEiDA_stats_FracOccup(leida_res,Conditions_tag,Paired_tests);

% Compute the dwell time and perform hypothesis tests
LEiDA_stats_DwellTime(leida_res,Conditions_tag,Paired_tests,TR);

% Generate and save the p-value and barplot plots for fractional occupancy
Plot_FracOccup(leida_res)

% Generate and save the p-value and barplot plots for dwell time
Plot_DwellTime(leida_res)

% Plot the centroids obtained using LEiDA and their overlap with Yeo nets
Plot_Centroid_Pyramid(leida_res,Conditions_tag,Parcellation,N_areas,CortexDirection)