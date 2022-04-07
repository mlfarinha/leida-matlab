function LEiDA_AnalysisCentroid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to analyse one PL state from the set of K PL states chosen
% according to the analysis from LEiDA_Start. This script will provide
% detailed figures with relevant information regarding the selected
% PL state. This function contains two sections: (A) user input parameters;
% and (B) code to plot the figures for the selected PL state. The user
% should only need to adapt section A.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Analysis plots for selected centroid:
%    - Plot PL state as vector with labels of parcellation
%    - Plot PL state as vector with areas ordered by contribution
%    - Plot boxplot of fractional occupancy values by condition
%    - Plot boxplot of dwell time values by condition
%    - Plot of summary information for the selected PL state
%
% Tutorial: README.md
% Version:  V1.0, April 2022
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: USER INPUT PARAMETERS

% Define K value, i.e., K returning the most significant differences between conditions:
SelectK = 17;
% Define the PL state to be studied (1 <= PL state <= SelectK):
Centroid = 15;

% Directory of the LEiDA toolbox folder:
LEiDA_directory = 'D:/LEiDA_Toolbox/';
% Name of the run to be used to create the folder to save the data:
run_name = 'ABIDE_dparsf_all_AAL116';
% Parcellation used to run LEiDA_Start script:
Parcellation = 'AAL116';

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

%% ANALYSE PL STATE OF INTEREST

% Close all open figures
close all;

% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)

% Directory with the results from LEiDA
leida_res = [LEiDA_directory 'LEiDA_Results_' run_name '/'];

% Check if the value for centroid is within the allowed range
if Centroid < 1 || Centroid > SelectK
    error(['Select a PL state within the range 1 to ' num2str(SelectK)])
end

% Create a directory to store results for defined value of K
if ~exist([leida_res 'K' num2str(SelectK) 'C' num2str(Centroid) '/'], 'dir')
    mkdir([leida_res 'K' num2str(SelectK) 'C' num2str(Centroid) '/']);
end
Vc_dir = [leida_res 'K' num2str(SelectK) 'C' num2str(Centroid) '/'];

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES FOR PL STATE ' num2str(Centroid) ' (K = ' num2str(SelectK) ') %%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% Plot PL state in vector format with label of each parcel
Plot_C_vector_labelled(leida_res,Vc_dir,SelectK,Centroid,Parcellation);

% Plot PL state in vector format with areas organised by contribution
Plot_C_vector_ordered(leida_res,Vc_dir,SelectK,Centroid,Parcellation);

% Plot boxplot of the fractional occupancy values for the selected PL state
Plot_C_boxplot_FO(leida_res,Vc_dir,SelectK,Centroid,Parcellation);

% Plot boxplot of the dwell time values for the selected PL state
Plot_C_boxplot_DT(leida_res,Vc_dir,SelectK,Centroid,Parcellation);

% Plot of summary information for the selected PL state
Plot_C_summary(leida_res,Vc_dir,SelectK,Centroid,Parcellation);