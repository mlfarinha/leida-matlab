function abide_get_data(data_dir,info_file,save_dir)
%
% Extract the time series for each parcel from the ABIDE dataset. The
% derivative rois_aal was extracted as .1D file. This script gets the
% time series data and stores it as a .mat file.
%
% INPUT:
% data_dir      directory where the .1D files with the ABIDE data are
%               stored
% info_file     path to file with the phenotypic data from each subject
% save_dir      directory to save the new data
%
% Author: Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% Input examples:
% data_dir = 'D:/LEiDA_Toolbox/Outputs/dparsf/nofilt_noglobal/rois_aal/';
% info_file = 'D:/LEiDA_Toolbox/Phenotypic_V1_0b.csv';
% save_dir = 'D:/LEiDA_Toolbox/ABIDE_dparsf_all_aal116/';

% Get number of files in folder
aux_data = [dir(fullfile([data_dir '*.mat'])); dir(fullfile([data_dir '*.1D'])); dir(fullfile([data_dir '*.txt']))];
num_subjs = numel(aux_data);

% Order the directory by name
[~,ind] = sort({aux_data.name});
data_info = aux_data(ind);

% Read the CSV with phenotypic data
info_data = readtable(info_file);
len_info = size(info_data,1);

n_ad = 0;
n_hc = 0;
tmax = 0;
for s = 1:num_subjs
    file = data_info(s).name;
    [~, baseFileName, ~] = fileparts(file);
    for i = 1:len_info
        pattern = num2str(info_data{i,2});
        if contains(baseFileName,pattern)
            disp(['File ' baseFileName ' parcellated using dparsf pipeline']);
            
            % Get data from .1D file
            data_file = importdata([data_dir file]);
            % Get data as N_areas*Tmax
            data = data_file.data';
            
            if tmax < size(data,2)
                tmax = size(data,2);
            end
        
            % If Diagnostic Group = autism (1 in table) and DSM-IV-TR Diagnostic Group =
            % Autism (1 in table) then append AD (Autism Disorder)
            if isequal(info_data{i,3},1) % && isequal(info_data{i,4},1)
                n_ad = n_ad + 1;
                disp('          - Participant with autism');
                save([save_dir baseFileName '_Autism'], 'data')
            % If Diagnostic Group = control (2 in table) and DSM-IV-TR Diagnostic Group =
            % control (0 in table) then append HC (Healthy Control)
            elseif isequal(info_data{i,3},2) % && isequal(info_data{i,4},0)
                n_hc = n_hc + 1;
                disp('          - Participant is healthy control');
                save([save_dir baseFileName '_Control'], 'data')
            else
                % disp('          - Participant with other autism spectrum disorders');
                disp('          - Group not specified');
                continue
            end         
        end
    end
end
disp(['Number of participants with tag Control: ' num2str(n_hc)]); % 476
disp(['Number of participants with tag Autism: ' num2str(n_ad)]); % 408
disp(['The maximum number of TRs across participants is: ' num2str(tmax)]); % 315

