function Plot_TransitionMatrix(data_dir,bestK)
%
% Plot the results from the hypothesis tests obtained from comparing the
% mean state-to-state transition probability between conditions
%
% INPUT:
% data_dir      directory with the results from running the hypothesis
%               tests on the state-to-state transition probabilities
% bestK         optimal K defined by the user
%
% OUTPUT:
% Fig1          plot of the estimated TPM for each condition
% Fig2          plot of the two-sided p-values obtained from the comparison
%               of the mean state-to-state transition probabilities for 
%               each pair of conditions
%
% Authors: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with results for the fractional occupancy  (output from LEiDA_stats_FracOccup.m)
file_TM = ['LEiDA_Stats_TransitionMatrix_K' num2str(bestK) '.mat'];

% Load required data:
load([data_dir file_TM],'cond','TMnorm','TM_pval2sided','Index_Conditions','effectsize','rangeK');

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Code below generates the figures relative to the fractional occupancy
disp('%%%%%%%%%%%%%%%%%%%% TRANSITION PROBABILITIES Figures %%%%%%%%%%%%%%%%%%%%')

%% FIGURE 1: PLOT OF THE MEAN TRANSITION MATRIX FOR EACH CONDITION

disp('Plotting mean transition probability matrix for each condition:')
TM_cond = cell(1,n_Cond);
for j = 1:n_Cond
    % Compute mean TPM omitting the NaN values
    TM_cond{j} = squeeze(mean(TMnorm(Index_Conditions == j,:,:),'omitnan'));
end

for i = 1:n_Cond
    Fig1 = figure('Position', get(0, 'Screensize'));
    mean_TM = TM_cond{i};
    imagesc(mean_TM, 'AlphaData', 0.8);
    for c_out = 1:bestK
        for c_in = 1:bestK
            caption = sprintf('%.2f', mean_TM(c_out,c_in));
            text(c_in,c_out, caption, 'FontSize', 15, 'Color', [1, 1, 1],'FontWeight','bold','HorizontalAlignment','Center');
        end
    end
    colormap('jet');
    set(colorbar,'YTick',[0 0.2 0.4 0.6 0.8 1], 'ylim', [0 1], 'FontSize', 12);
    ticklabels = cell(1,bestK);
    for p = 1:bestK
        ticklabels{p} = p;
    end
    xticks(1:1:bestK); xticklabels(ticklabels);
    yticks(1:1:bestK); yticklabels(ticklabels);
    set(gca,'FontSize',14,'FontWeight','bold')
    axis square
    ylabel('From FC State')
    xlabel('To FC State')
    
    saveas(Fig1, fullfile(data_dir, ['OptK' num2str(bestK) '_TPM_' cond{i} '.png']),'png');
    saveas(Fig1, fullfile(data_dir, ['OptK' num2str(bestK) '_TPM_' cond{i} '.fig']),'fig');
    disp(['- Mean transition probability matrix of condition ' cond{i} ' successfully saved as OptK' num2str(bestK) '_TPM_' cond{i}]);
end

close all;

%% FIGURE 2: PLOT OF DIFFERENCES BETWEEN TRANSITION PROBABILITIES ACROSS CONDITIONS

% Create a map to place the figures
pos = power(n_Cond-1,2);
index_fig = reshape(1:pos, n_Cond-1, n_Cond-1).';
subplot_map = ones(n_Cond-1);
subplot_map = triu(subplot_map).';

% Define indices for conditions in subplot
subplot_indices = find(subplot_map);

% Possible pairs of conditions comparisons
condRow = zeros(1,n_Cond*(n_Cond-1)/2);
condCol = zeros(1,n_Cond*(n_Cond-1)/2);
cond_pair = 1;
% disp('Possible pairs of condition comparisons:')
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond
        condRow(cond_pair) = cond1;
        condCol(cond_pair) = cond2;
        % disp([num2str(cond{cond1}) ' - ' num2str(cond{cond2})])
        cond_pair = cond_pair + 1;
    end
end

Fig2 = figure('Position', get(0, 'Screensize'));
disp('Plotting two-sided p-values from hypothesis tests:')
for s_ind = 1:length(subplot_indices)
    
    % disp(['- ' cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}])
    
    subplot_ind = subplot_indices(s_ind);
    
    % Subplot adds plots following the rows (not the columns)
    subplot(size(subplot_map,1),size(subplot_map,2),index_fig(subplot_ind))
    
    % Something in the background
    I = magic(bestK);
    % Generate where each text will go
    [X, Y] = meshgrid(1:bestK,1:bestK);
    % Display the background
    imagesc(I);
    hold on;
    % Set the background color to white
    set(gcf,'Colormap',[1, 1, 1]);
    % Insert the labels
    for c_out = 1:bestK
        for c_in = 1:bestK
            if TM_pval2sided(s_ind,c_out,c_in) < (0.05/bestK) && TM_pval2sided(s_ind,c_out,c_in) > (0.05/sum(rangeK))
                text(X(c_in,c_out),Y(c_out,c_in)+.22,'*','FontSize',18,'HorizontalAlignment','center','FontWeight','bold','Color','g');
            end
            if TM_pval2sided(s_ind,c_out,c_in) <= (0.05/sum(rangeK))
                text(X(c_in,c_out),Y(c_out,c_in)+.22,'*','FontSize', 18,'HorizontalAlignment','center','FontWeight','bold','Color','b');
            end
        end
    end
    % Calculte the grid lines
    grid = 0.5:1:bestK+0.5;
    grid1 = [grid; grid];
    grid2 = repmat([0.5;bestK+0.5],1,length(grid));
    % Plot the grid lines
    plot(grid1,grid2,'k')
    plot(grid2,grid1,'k')
    ticklabels = cell(1,bestK);
    for p = 1:bestK
        ticklabels{p} = p;
    end
    xticks(1:1:bestK); xticklabels(ticklabels);
    yticks(1:1:bestK); yticklabels(ticklabels);
    set(gca,'FontSize',9)
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}],'interpreter','none')
    ylabel('From FC State','FontSize',11)
    xlabel('To FC State','FontSize',11)
    axis square
    box off
    
    sgtitle([{'Two-sided {\itp}-value'},{'State-to-state transition probability'}],'Fontsize',14,'FontWeight','bold') 
end
saveas(Fig2, fullfile(data_dir, ['OptK' num2str(bestK) '_TransProbs_pvalues.png']),'png');
saveas(Fig2, fullfile(data_dir, ['OptK' num2str(bestK) '_TransProbs_pvalues.fig']),'fig');
disp('- Plot successfully saved as TransProbs_pvalues');

close all;