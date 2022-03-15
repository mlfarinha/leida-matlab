function LEiDA_OptimalK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            LEADING EIGENVECTOR DYNAMICS ANALYSIS - OPTIMAL K            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function runs the LEiDA analysis for the optimal number of FC states
% selected based on the outputs from LEiDA_Master.m. This function contains
% two sections: (A) user input parameters; and (B) code to run LEiDA for
% optimal K. The user should only need to adapt section A.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Run Leading Eigenvector Dynamics Analysis for optimal K:
%    - Plot several plots for the selected optimal K
%    - Analysis of the transition probability matrix values
%
% Tutorial: README.md
% Version:  V1.0, March 2022
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: USER INPUT PARAMETERS

% Define optimal K value, i.e., optimal number of FC states:
OptimalK = 10;

% Directory of the LEiDA toolbox folder:
LEiDA_directory = 'D:\LEiDA_Toolbox\';
% Parcellation applied to the imaging data:
Parcellation = 'AAL116';
% Number of brain areas to consider for analysis;
N_areas = 90;
% Experimental paradigm (0: subjects in different conditions are not the
% same; 1: subjects are the same across conditions)
Paired_tests = 0;

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

%% RUN ANALYSIS FOR OPTIMAL K

% Close all open figures
close all;

% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)

% Directory with the results from LEiDA
leida_res = [LEiDA_directory '\LEiDA_Results\'];

% Create a directory to store results for defined value of FC states
mkdir([leida_res num2str(OptimalK) 'FCstates\']);
OptK_dir = [leida_res num2str(OptimalK) 'FCstates\'];

% Summary plot of the optimal K
Plot_OptimalK(leida_res,OptK_dir,OptimalK,Parcellation,N_areas);

% Compute the transition probability matrix and perform hypothesis tests
LEiDA_stats_TransitionMatrix(leida_res,OptimalK,Paired_tests);

% Plot the mean transition probability matrix and the differences in
% state-to-state transition probabilities between conditions
Plot_TransitionMatrix(leida_res,OptK_dir,OptimalK);

