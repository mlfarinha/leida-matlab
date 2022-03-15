function Parcellate


%% USER INPUT PARAMETERS

% Directory of the folder with the neuroimaging data:
fmri_dir = 'D:\IST\tese\res\noNR_data\';
% Directory to save the parcellated data:
save_dir = 'D:\LEiDA_Toolbox\fMRIdata_parcellated_AAL116\';
% Define the parcellation you would like to apply to the data:
% Available parcellations: AAL116, dbs80, glasser360, glasser_378, Yeo7
Parcellation = 'AAL116';
% Define the number of areas in which the brain will be parcellated:
N_areas = 90;

%% RUN TO PARCELLATE THE NEUROIMAGING DATA

% Get number of files in folder
aux_data = dir([fmri_dir '*.nii']);
num_subjs = numel(aux_data);

% Obtain parcellation atlas from ParcelsMNI2mm (Dr. Joana Cabral)
V_Parcel = struct2array(load('ParcelsMNI2mm',['V_' Parcellation])); 
sz = size(V_Parcel); % size of the parcellation

for s = 1:num_subjs
    
    file = aux_data(s).name;
    [~, baseFileName, ~] = fileparts(file);
    
    if size(file,1)
        
        disp(['Parcellating data file ' baseFileName ' using parcellation ' Parcellation]);
        
        % Read the nii file
        fMRI_MNI = niftiread([fmri_dir file]);
        T = size(fMRI_MNI,4); % number of volumes
        
        % Check if nii files correspond to MNI2mm
        if size(fMRI_MNI,1) ~= sz(1) || size(fMRI_MNI,2) ~= sz(2) || size(fMRI_MNI,3) ~= sz(3)
            % Files will be resized in order to be accomodated to the MNI2mm template
            fMRI_MNI2mm = zeros(sz(1), sz(2), sz(3), T);
            for t=1:T
                fMRI_MNI2mm(:,:,:,t) = imresize3(fMRI_MNI(:,:,:,t),sz);
            end
        end
                
        fMRI_parcel = zeros(N_areas,T);
        
        for n = 1:N_areas
            ind_voxels = find(V_Parcel==n);
            
            for v = 1:numel(ind_voxels)
                [I1,I2,I3] = ind2sub(sz,ind_voxels(v));
                if ~isnan(fMRI_MNI2mm(I1,I2,I3,1))
                    fMRI_parcel(n,:) = fMRI_parcel(n,:) + squeeze(fMRI_MNI2mm(I1,I2,I3,:))';
                end
            end

            fMRI_parcel(n,:) = fMRI_parcel(n,:)/numel(ind_voxels);
            fMRI_parcel(n,:) = detrend(fMRI_parcel(n,:) - mean(fMRI_parcel(n,:)));
            fMRI_parcel(n,:) = fMRI_parcel(n,:)/std(fMRI_parcel(n,:));
        end
        
        save([save_dir baseFileName '_' Parcellation], 'fMRI_parcel')
    end
end  