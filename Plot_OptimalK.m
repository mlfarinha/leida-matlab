function Plot_OptimalK(data_dir,save_dir,bestK,parcellation,n_areas)
%
% Plot the results from the hypothesis tests obtained from comparing the
% mean fractional occupancy between conditions
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% bestK         optimal K defined by the user
% parcellation  parcellation template used to segment the brain
% n_areas       number of areas from the parcellation to analyse
%
% OUTPUT:
% Fig1          Set of FC states rendered in 3D glass brain
% Fig2          Set of FC states in matrix format
% Fig3          Set of FC states in vector format with labels of
%               parcellation
% Fig4          Set of FC states in vector format with numbered parcels
% Fig5          FC states in vector format organised by contribution
% Fig6          Boxplots of Fractional Occupancy values of all FC states 
% Fig7          Boxplot of Fractional Occupancy values for each FC state
% Fig8          Boxplots of Dwell Time values of all FC states
% Fig9          Boxplot of Dwell Time values for each FC state
% Fig10         Plot of vector, 3D glass brain, matrix and boxplots of
%               each FC state
% Fig11         Plot of 3D glass brain, vector and boxplots of all
%               FC states
%
% Authors: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_DwellTime.m)
file_LT = 'LEiDA_Stats_DwellTime.mat';

% Load required data:
load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
load([data_dir file_P], 'cond', 'P', 'P_pval2sided', 'Index_Conditions');
load([data_dir file_LT], 'LT', 'LT_pval2sided');

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Color code from paper Yeo et al., 2011
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;
Yeo_names = {'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);

% LEiDA networks colored according to closest RSN
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation]));
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Matrix of dimension bestK*90 where each row represents one FC state
V = Kmeans_results{rangeK == bestK}.C;
% Scale each cluster centroid by its maximum value and transpose the matrix
V = V'./max(abs(V'));
clear Kmeans_results

% Get labels from parcellation used
switch parcellation
    case 'AAL116'
        load AAL_labels.mat label90
        labels_parcel = label90([1:n_areas],:);
        clear label90
end

% Reorder the position of each parcel
Order = [1:2:n_areas n_areas:-2:2];

% Code below generates the figures relative to the fractional occupancy
disp('%%%%%%%%%%%%%%%%%%%% Optimal K Figures %%%%%%%%%%%%%%%%%%%%')

%% FIGURE 1: FC STATE RENDERED IN 3D GLASS BRAIN

disp('Plotting FC states in 3D glass brain:')
Fig1 = figure('Position', get(0, 'Screensize'));
for c = 1:bestK
    
    [~, net] = max(cc_V_yeo7(rangeK == bestK,c,:));
    
    % First Plot view from top
    subplot_tight(2,bestK,c,0.01)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Color all areas with positive V with transparency proportional to
    % contribution
    n_pos = find(V(:,c) > 0);
    if numel(n_pos) > 0
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else % Plot of the global mode with "gray/transparent" color
        region_Vol = Volume > 0 & Volume <= 90; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        title({['FC state ' num2str(c), newline, Yeo_names{net}]}, 'Fontsize', 12) %,Yeo_names{c}})
    elseif numel(n_pos) == 0
        title({['FC state ' num2str(c), newline, 'Global Mode']}, 'Fontsize', 12) %,Yeo_names{c}})
    else
        title({['FC state ' num2str(c), newline, '']}, 'Fontsize', 12) %,Yeo_names{c}})
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %  Same but view from the side
    subplot_tight(2,bestK,c+bestK,0.01)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    n_pos = find(V(:,c) > 0);
    if n_pos
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else
        region_Vol = Volume > 0 & Volume <= 90; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
end
saveas(Fig1, fullfile(save_dir, ['OptK' num2str(bestK) '_3DbrainRendering.png']),'png');
saveas(Fig1, fullfile(save_dir, ['OptK' num2str(bestK) '_3DbrainRendering.fig']),'fig');
disp(['- Plot successfully saved as OptK' num2str(bestK) '_3DbrainRendering']);

close all;
    

%% FIGURE 2: FC STATE IN MATRIX FORMAT

disp('Plotting FC states in matrix format:')
Fig2 = figure('Position', get(0, 'Screensize'));  
for c = 1:bestK
    subplot_tight(1,bestK,c,0.02)
    colormap(jet)
    imagesc(V(Order,c)*V(Order,c)',[-1 1])
    axis square
    title(['V_{C_{' num2str(c) '}}' '.V_{C_{' num2str(c) '}}' '^T'], 'Fontsize', 12)
end
saveas(Fig2, fullfile(save_dir, ['OptK' num2str(bestK) '_Matrix.png']),'png');
saveas(Fig2, fullfile(save_dir, ['OptK' num2str(bestK) '_Matrix.fig']),'fig');
disp(['- Matrices successfully saved as OptK' num2str(bestK) '_Matrix']);

close all;


%% FIGURE 3: FC STATE IN VECTOR FORMAT WITH LABELS OF PARCELLATION

disp('Plotting barplots of FC states in vector format with labels of parcellation:')
Fig3 = figure('Position', get(0, 'Screensize'));   
for c = 1:bestK
    subplot(1,bestK,c)
    Vo = V(Order,c);
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 n_areas+1])
    xlim([-1 1])
    set(gca,'YTick',1:n_areas)
    set(gca,'YTickLabel',labels_parcel(Order,:),'Fontsize',6)
    ax = gca;
    ax.XAxis.FontSize = 8;
    grid on
    title(['FC state ' num2str(c)],'Fontsize',10)
end
saveas(Fig3, fullfile(save_dir, ['OptK' num2str(bestK) '_VectorLabels' parcellation '.png']),'png');
saveas(Fig3, fullfile(save_dir, ['OptK' num2str(bestK) '_VectorLabels' parcellation '.fig']),'fig');
disp(['- Barplots successfully saved as OptK' num2str(bestK) '_VectorLabels' parcellation]);

close all;


%% FIGURE 4: FC STATES IN VECTOR FORMAT WITH NUMBERED PARCELS

disp('Plotting barplots of FC states in vector format with numbered parcels:')
Fig4 = figure('Position', get(0, 'Screensize'));   
for c = 1:bestK
    subplot(1,bestK,c)
    Vo = V(Order,c);
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 n_areas+1])
    xlim([-1 1])
    set(gca,'YTick',[10:10:n_areas], 'Fontsize',12)
    ax = gca;
    ax.XAxis.FontSize = 12;
    grid on
    title(['FC state ' num2str(c)],'Fontsize',12)
end
saveas(Fig4, fullfile(save_dir, ['OptK' num2str(bestK) '_VectorNumber.png']),'png');
saveas(Fig4, fullfile(save_dir, ['OptK' num2str(bestK) '_VectorNumber.fig']),'fig');
disp(['- Barplots successfully saved as OptK' num2str(bestK) '_VectorNumber']);

close all;


%% FIGURE 5: INDIVIDUAL FC STATES IN VECTOR FORMAT ORGANISED BY CONTRIBUTION

disp('Plotting barplot of each FC state individually:')  
for c = 1:bestK
    [Vo, asc_ind] = sort(V(:,c));
    numMarkers = size(Vo,1);
    if c == 1
        markerColors = jet(numMarkers*2);
    else
        markerColors = jet(numMarkers);
    end
    miny = min(Vo);
    maxy = max(Vo);
    colorMapRows = round(rescale((Vo - miny) / (maxy - miny), 1, numMarkers));
    Fig5 = figure('Position', get(0, 'Screensize'));
    subplot(1,3,2)
    for k = 1:length(Vo)
        thisMarkerColor = markerColors(colorMapRows(k), :);
        barh(k,Vo(k),'FaceColor',thisMarkerColor,'EdgeColor','none','Barwidth',.5)
        hold on;
        ylim([0 n_areas+1])
        xlim([-1 1])
        set(gca,'XTick',[-1 0 1])
        set(gca,'YTick',1:n_areas)
        set(gca,'YTickLabel',labels_parcel(asc_ind,:),'Fontsize',6)
        ax = gca;
        ax.XAxis.FontSize = 10;
        grid on;
        title(['FC state ' num2str(c)],'Fontsize',12)
    end
    saveas(Fig5, fullfile(save_dir, ['OptK' num2str(bestK) '_FCstate' num2str(c) '.png']),'png');
    saveas(Fig5, fullfile(save_dir, ['OptK' num2str(bestK) '_FCstate' num2str(c) '.fig']),'fig');
    disp(['- Barplot of FC state ' num2str(c) ' successfully saved as OptK' num2str(bestK) '_FCstate' num2str(c)]);
end

close all;


%% FIGURE 6: BOXPLOT OF THE FRACTIONAL OCCUPANCY VALUES OF ALL FC STATES

disp('Plotting boxplots of the Fractional Occupancy values:')
Fig6 = figure('Position', get(0, 'Screensize'));   
for c = 1:bestK
    subplot(4,bestK,[c+bestK c+2*bestK])
    [~, net] = max(cc_V_yeo7(rangeK == bestK,c,:));
    P_cond = cell(1,n_Cond);
    g = cell(1,n_Cond);
    for j = 1:n_Cond
        P_cond{j} = P(Index_Conditions == j,rangeK == bestK,c);
        g{j} = repmat(cond(j), length(P(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    P_data = vertcat(P_cond{:});
    g_data = vertcat(g{:});
    
    bx = boxplot(P_data, g_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end   
    
    hold on
    set(gca,'XTickLabel',cond)
    xtickangle(30)
    if c == 1
        ylabel('Fractional Occupancy (%)', 'Fontsize', 10);
    end
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(P_data);
    Y_LOCATIONS = Max_Y + 0.05*(abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    % end
    
    % Blue asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    % end
    
    ylim([0 max(Y_LOCATIONS)+0.15])
    
    hold off
    box off
end
saveas(Fig6, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotFracOccup.png']),'png');
saveas(Fig6, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotFracOccup.fig']),'fig');
disp(['- Boxplots successfully saved as OptK' num2str(bestK) '_BoxPlotFracOccup']);

close all;


%% FIGURE 7: BOXPLOT OF THE FRACTIONAL OCCUPANCY VALUES OF EACH FC STATE

disp('Plotting boxplot of the Fractional Occupancy values for each FC state:')  
for c = 1:bestK
    [~, net] = max(cc_V_yeo7(rangeK == bestK,c,:));
    P_cond = cell(1,n_Cond);
    g = cell(1,n_Cond);
    for j = 1:n_Cond
        P_cond{j} = P(Index_Conditions == j,rangeK == bestK,c);
        g{j} = repmat(cond(j), length(P(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    P_data = vertcat(P_cond{:});
    g_data = vertcat(g{:});
    
    Fig7 = figure;
    bx = boxplot(P_data, g_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end   
    
    hold on
    set(gca,'XTickLabel',cond,'Fontsize',12)
    xtickangle(30)
    ylabel('Fractional Occupancy (%)','Fontsize',12);
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(P_data);
    Y_LOCATIONS = Max_Y + 0.05*(abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)

    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    % Blue asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    ylim([0 max(Y_LOCATIONS)+0.15])
    
    hold off
    box off
    
    saveas(Fig7, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotFracOccup_FCstate' num2str(c) '.png']),'png');
    saveas(Fig7, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotFracOccup_FCstate' num2str(c) '.fig']),'fig');
    disp(['- Boxplot successfully saved as OptK' num2str(bestK) '_BoxPlotFracOccup_FCstate' num2str(c)]);
end

close all;


%% FIGURE 8: BOXPLOT OF THE DWELL TIME VALUES OF ALL FC STATES

disp('Plotting boxplots of the Dwell Time values:')
Fig8 = figure('Position', get(0, 'Screensize'));      
for c = 1:bestK
    subplot(4,bestK,[c+bestK c+2*bestK])
    [~, net] = max(cc_V_yeo7(rangeK == bestK,c,:));
    LT_cond = cell(1,n_Cond);
    r = cell(1,n_Cond);
    for j = 1:n_Cond
        LT_cond{j} = LT(Index_Conditions == j,rangeK == bestK,c);
        r{j} = repmat(cond(j), length(LT(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    LT_data = vertcat(LT_cond{:});
    r_data = vertcat(r{:});
    
    bx = boxplot(LT_data, r_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    z = get(get(gca,'children'),'children');   % Get the handles of all the objects
    y = get(z,'tag');   % List the names of all the objects
    med_col = z(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = z(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end
    
    hold on
    set(gca,'XTickLabel',cond)
    xtickangle(30)
    if c == 1
        ylabel('Dwell Time (s)','Fontsize',10);
    end
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(LT_data);
    Y_LOCATIONS = Max_Y + (abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    % end
    
    % Blue asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    % end
    
    ylim([0 max(Y_LOCATIONS)*1.15])
    
    hold off
    box off
end
saveas(Fig8, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotDwellTime.png']),'png');
saveas(Fig8, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotDwellTime.fig']),'fig');
disp(['- Boxplots successfully saved as OptK' num2str(bestK) '_BoxPlotDwellTime']);

close all;


%% FIGURE 9: BOXPLOT OF THE DWELL TIME VALUES OF EACH FC STATE

disp('Plotting boxplot of the Dwell Time values for each FC state:')  
for c = 1:bestK
    [~, net] = max(cc_V_yeo7(rangeK == bestK,c,:));
    LT_cond = cell(1,n_Cond);
    r = cell(1,n_Cond);
    for j = 1:n_Cond
        LT_cond{j} = LT(Index_Conditions == j,rangeK == bestK,c);
        r{j} = repmat(cond(j), length(LT(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    LT_data = vertcat(LT_cond{:});
    r_data = vertcat(r{:});
    
    Fig9 = figure;
    bx = boxplot(LT_data, r_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end
    
    hold on
    set(gca,'XTickLabel',cond,'Fontsize',12)
    xtickangle(30)
    ylabel('Dwell Time (s)','Fontsize',12);
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(LT_data);
    Y_LOCATIONS = Max_Y + (abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)

    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    % Blue asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    ylim([0 max(Y_LOCATIONS)*1.15])
    
    hold off
    box off
    
    saveas(Fig9, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotDwellTime_FCstate' num2str(c) '.png']),'png');
    saveas(Fig9, fullfile(save_dir, ['OptK' num2str(bestK) '_BoxPlotDwellTime_FCstate' num2str(c) '.fig']),'fig');
    disp(['- Boxplot successfully saved as OptK' num2str(bestK) '_BoxPlotDwellTime_FCstate' num2str(c)]);
end

close all;


%% FIGURE 10: PLOT OF VECTOR, 3D GLASS BRAIN, MATRIX AND BOXPLOTS FOR EACH FC STATE

disp('Plotting vector, 3D glass brain and matrix and boxplots for each FC state:')
for c = 1:bestK
    Fig10 = figure('Position', get(0, 'Screensize'));  
    
    % Plot of vector with parcels organized by contribution
    [Vo, asc_ind] = sort(V(:,c));
    numMarkers = size(Vo,1);
    if c == 1
        markerColors = jet(numMarkers*2);
    else
        markerColors = jet(numMarkers);
    end
    miny = min(Vo);
    maxy = max(Vo);
    colorMapRows = round(rescale((Vo - miny) / (maxy - miny), 1, numMarkers));
    subplot(2,4,[1 5])
    for k = 1:length(Vo)
        thisMarkerColor = markerColors(colorMapRows(k), :);
        barh(k,Vo(k),'FaceColor',thisMarkerColor,'EdgeColor','none','Barwidth',.5)
        hold on;
        ylim([0 n_areas+1])
        xlim([-1 1])
        set(gca,'XTick',[-1 0 1])
        set(gca,'XTickLabel',[-1 0 1])
        set(gca,'YTick',1:n_areas)
        set(gca,'YTickLabel',labels_parcel(asc_ind,:),'Fontsize',6)
        ax = gca;
        ax.XAxis.FontSize = 10;
        grid on;
        title({'fMRI phase','projection into V_{1}'},'Fontsize',10)
    end
    
    % PLOT OF 3D GLASS BRAIN
    [~, net] = max(cc_V_yeo7(rangeK == bestK,c,:));
    subplot(2,4,2)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Color all areas with positive V with transparency proportional to
    % contribution
    n_pos = find(V(:,c) > 0);
    if numel(n_pos) > 0
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else % Plot of the global mode with "gray/transparent" color
        region_Vol = Volume > 0 & Volume <= 90; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    title(['FC state ' num2str(c)], 'Fontsize', 10)    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    
    % PLOT OF MATRIX FORMAT
    subplot(2,4,3)
    colormap(jet)
    imagesc(V(Order,c)*V(Order,c)',[-1 1])
    axis square
    title(['V_{C_{' num2str(c) '}}' '.V_{C_{' num2str(c) '}}' '^T'], 'Fontsize', 10)
    
    
    % PLOT OF FRACTIONAL OCCUPANCY VALUES
    subplot(2,4,6)
    P_cond = cell(1,n_Cond);
    g = cell(1,n_Cond);
    for j = 1:n_Cond
        P_cond{j} = P(Index_Conditions == j,rangeK == bestK,c);
        g{j} = repmat(cond(j), length(P(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    P_data = vertcat(P_cond{:});
    g_data = vertcat(g{:});
    
    bx = boxplot(P_data, g_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end   
    
    hold on
    set(gca,'XTickLabel',cond,'Fontsize',10)
    xtickangle(30)
    ylabel('Fractional Occupancy (%)','Fontsize',10);
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(P_data);
    Y_LOCATIONS = Max_Y + 0.05*(abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)

    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    % Blue asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    ylim([0 max(Y_LOCATIONS)+0.15])
    
    hold off
    box off
    
    % PLOT OF DWELL TIME VALUES
    subplot(2,4,7)
    LT_cond = cell(1,n_Cond);
    r = cell(1,n_Cond);
    for j = 1:n_Cond
        LT_cond{j} = LT(Index_Conditions == j,rangeK == bestK,c);
        r{j} = repmat(cond(j), length(LT(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    LT_data = vertcat(LT_cond{:});
    r_data = vertcat(r{:});
    
    bx = boxplot(LT_data, r_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end
    
    hold on
    set(gca,'XTickLabel',cond,'Fontsize',10)
    xtickangle(30)
    ylabel('Dwell Time (s)','Fontsize',10);
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(LT_data);
    Y_LOCATIONS = Max_Y + (abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)

    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    % Blue asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    if asterisks
        text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    end
    
    ylim([0 max(Y_LOCATIONS)*1.15])
    
    hold off
    box off
    
    saveas(Fig10, fullfile(save_dir, ['OptK' num2str(bestK) '_Summary_FCstate' num2str(c) '.png']),'png');
    saveas(Fig10, fullfile(save_dir, ['OptK' num2str(bestK) '_Summary_FCstate' num2str(c) '.fig']),'fig');
    disp(['- Summary plot of FC state ' num2str(c) ' successfully saved as OptK' num2str(bestK) '_Summary_FCstate' num2str(c)]);
end

close all;


%% FIGURE 11: PLOT OF 3D GLASS BRAIN, VECTOR AND BOXPLOTS

disp('Plotting 3D glass brain, vector and boxplots:')
Fig11 = figure('Position', get(0, 'Screensize')); 
colormap(jet)
% set(gcf, 'Position',  [100, 100, 1000, 500])
for c = 1:bestK
    
    [~, net] = max(cc_V_yeo7(rangeK == bestK,c,:));
    
    % First Plot view from top
    subplot(9,bestK,c)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Color all areas with positive V with transparency proportional to
    % contribution
    n_pos = find(V(:,c) > 0);
    if numel(n_pos) > 0
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else % Plot of the global mode with "gray/transparent" color
        region_Vol = Volume > 0 & Volume <= 90; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        title({['FC state ' num2str(c), newline, Yeo_names{net}]} ) %,Yeo_names{c}})
    elseif numel(n_pos) == 0
        title({['FC state ' num2str(c), newline, 'Global Mode']} ) %,Yeo_names{c}})
    else
        title({['FC state ' num2str(c), newline, '']} ) %,Yeo_names{c}})
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %  Same but view from the side
    subplot(9,bestK,c+bestK)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    n_pos = find(V(:,c) > 0);
    if n_pos
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else
        region_Vol = Volume > 0 & Volume <= 90; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    % PLOT THE VECTOR OF EACH FC STATE
    subplot(9,bestK,[c+2*bestK c+3*bestK c+4*bestK])
    Vo = V(Order,c);
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 n_areas+1])
    xlim([-1 1])
    set(gca,'YTick',[10:20:n_areas], 'Fontsize',8)
    ax = gca;
    ax.XAxis.FontSize = 8;
    grid on
    
    % BOXPLOT FRACTIONAL OCCUPANCY
    subplot(9,bestK,[c+5*bestK c+6*bestK])
    
    P_cond = cell(1,n_Cond);
    g = cell(1,n_Cond);
    for j = 1:n_Cond
        P_cond{j} = P(Index_Conditions == j,rangeK == bestK,c);
        g{j} = repmat(cond(j), length(P(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    P_data = vertcat(P_cond{:});
    g_data = vertcat(g{:});
    
    bx = boxplot(P_data, g_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end   
    
    hold on
    set(gca,'XTickLabel',cond, 'Fontsize',6)
    xtickangle(30)
    if c == 1
        ylabel('Fractional Occupancy (%)');
    end
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(P_data);
    Y_LOCATIONS = Max_Y + 0.05*(abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    % end
    
    % Blue asterisks
    asterisks = find(P_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    % end
    
    ylim([0 max(Y_LOCATIONS)+0.15])
    
    hold off
    box off
    
    % BOXPLOT DWELL TIME
    subplot(9,bestK,[c+7*bestK c+8*bestK])
    
    LT_cond = cell(1,n_Cond);
    r = cell(1,n_Cond);
    for j = 1:n_Cond
        LT_cond{j} = LT(Index_Conditions == j,rangeK == bestK,c);
        r{j} = repmat(cond(j), length(LT(Index_Conditions == j,rangeK == bestK,c)),1);
    end
    LT_data = vertcat(LT_cond{:});
    r_data = vertcat(r{:});
    
    bx = boxplot(LT_data, r_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == bestK,c,net) < 0.05/bestK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end
    
    hold on
    set(gca,'XTickLabel',cond, 'Fontsize',6)
    xtickangle(30)
    if c == 1
        ylabel('Dwell Time (s)');
    end
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(LT_data);
    Y_LOCATIONS = Max_Y + (abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/bestK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
    % end
    
    % Blue asterisks
    asterisks = find(LT_pval2sided(:,rangeK == bestK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
     % if asterisks
         % text(mean(X_locations(asterisks,:),2) - 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
             % [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == bestK,c),'%10.1e')],'Fontsize',7)
     % end
    
    ylim([0 max(Y_LOCATIONS)*1.15])
    
    hold off
    box off
end
saveas(Fig11, fullfile(save_dir, ['OptK' num2str(bestK) '_Summary' num2str(c) '.png']),'png');
saveas(Fig11, fullfile(save_dir, ['OptK' num2str(bestK) '_Summary' num2str(c) '.fig']),'fig');
disp(['- Summary plot of all FC states successfully saved as OptK' num2str(bestK) '_Summary']);

close all;   
