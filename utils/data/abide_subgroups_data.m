function abide_subgroups_data(data_dir,info_file,save_dir)
%
% Extract the time series for each parcel from the ABIDE dataset. The
% derivative rois_aal was extracted as .1D file. This script gets the
% time series data, stores it as a .mat file and adds a tag to the filename
% of each participant to determine the group.
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
% save_dir = 'D:/LEiDA_Toolbox/ABIDE_dparsf_subconds_aal116/';

% Get number of files in folder
aux_data = [dir(fullfile([data_dir '*.mat'])); dir(fullfile([data_dir '*.1D'])); dir(fullfile([data_dir '*.txt']))];
num_subjs = numel(aux_data);

% Order the directory by name
[~,ind] = sort({aux_data.name});
data_info = aux_data(ind);

% Read the CSV with phenotypic data
info_data = readtable(info_file);
len_info = size(info_data,1);

% Count number of subjects from each group
n_ad = 0;
n_hc = 0;
n_asp = 0;
n_pdd = 0;
n_asp_or_pdd = 0;
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
            
            % Column label: DSM_IV_TR (column 4)
            % (1) Control == 0
            if isequal(info_data{i,4},0)
                n_hc = n_hc + 1;
                disp('          - Participant -> CONTROL');
                save([save_dir baseFileName '_CONT'], 'data')
            % (2) Autism == 1
            elseif isequal(info_data{i,4},1)
                n_ad = n_ad + 1;
                disp('          - Participant -> AUTISM');
                save([save_dir baseFileName '_AUT'], 'data')
            % (3) Aspergers = 2
            elseif isequal(info_data{i,4},2)
                n_asp = n_asp + 1;
                disp('          - Participant -> ASPERGERS');
                save([save_dir baseFileName '_ASP'], 'data')
            % (4) PDD-NOS = 3
            elseif isequal(info_data{i,4},3)
                n_pdd = n_pdd + 1;
                disp('          - Participant -> PDD-NOS');
                save([save_dir baseFileName '_PDD-NOS'], 'data')
            % (5) Aspergers or PDD-NOS = 4
            elseif isequal(info_data{i,4},4)
                n_asp_or_pdd = n_asp_or_pdd + 1;
                disp('          - Participant -> ASPERGERS OR PDD-NOS');
                save([save_dir baseFileName '_ASP-PDD'], 'data')
            else
                disp('          - NOT SPECIFIED');
                continue
            end         
        end
    end
end
disp(['Number of participants with tag Control: ' num2str(n_hc)]); % 460
disp(['Number of participants with tag Autism: ' num2str(n_ad)]); % 251
disp(['Number of participants with tag Aspergers: ' num2str(n_asp)]); % 72
disp(['Number of participants with tag PDD-NOS: ' num2str(n_pdd)]); % 32
disp(['Number of participants with tag Aspergers or PDD-NOS: ' num2str(n_asp_or_pdd)]); % 6
disp(['The maximum number of TRs across participants is: ' num2str(tmax)]); % 315