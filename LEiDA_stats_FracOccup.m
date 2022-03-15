function LEiDA_stats_FracOccup(data_dir,cond,pair)
%
% For each clustering solution compute the fractional occupancy of each FC
% state for all subjects. Perform hypothesis tests to check for differences
% in the probability of occurrence of FC states between conditions.
%
% INPUT:
% data_dir      directory where the LEiDA results are saved
% cond          tags of each condition considered in the experiment
% pair          0, different subjects in each condition (default); 1, same
%               subjects across conditions
%
% OUTPUT:
% P             fractional occupancy of each FC state in each clustering
%               solution for all subjects
% P_pval2sided  two-sided p-value from testing whether the mean fractional
%               occupancy of a given FC state differs between conditions
% effectsize    Hedge's effect size
% levene_pval   p-value obtained from computing the Levene's test on the
%               fractional occupancy values
%
% Author: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% Default number of permutations
n_permutations = 5000;
% Default number of bootstrap samples within each permutation sample
n_bootstraps = 500;

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with K-means results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
load([data_dir file_V1],'Time_sessions','Data_info');
load([data_dir file_cluster],'Kmeans_results', 'rangeK');

N_scans = max(Time_sessions);

% Number of conditions of the experiment
n_Cond = size(cond,2);

Index_Conditions = [];
for s = 1:length(Data_info)
    FileName = Data_info(s).name;
    for c = 1:n_Cond
        if contains(FileName,string(cond(c)))
            Index_Conditions = cat(2, Index_Conditions, c);
        else
            continue
        end
    end
end

% Code below computes the fractional occupancy and runs hypothesis tests
disp('%%%%%%%%%%%%%%%%%%%%%%% Fractional Occupancy %%%%%%%%%%%%%%%%%%%%%%%')

%% CALCULATE THE FRACTIONAL OCCUPANCY OF EACH STATE FOR K-MEANS CLUSTERING SOLUTIONS

% For every fMRI subject and condition calculate probability of each state c.
P = zeros(N_scans,length(rangeK),rangeK(end));

for k = 1:length(rangeK)
    for s = 1:N_scans

        T = Time_sessions == s;
        % Vector of the cluster to which each observation belongs to for subject s
        Ctime = Kmeans_results{k}.IDX(T);
        
        for c = 1:rangeK(k)
            % Probability of occurrence
            P(s,k,c) = mean(Ctime == c);
        end
    end
end

clear Kmeans_results

%% CHECK HOMOGENEITY OF VARIANCES USING LEVENE'S TEST

if pair == 0
    
    disp('Running Levene''s test of equality of variances:')
    % matrix to store p-value of Levene's test of equality of variances
    levene_pval = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));

    for k = 1:length(rangeK)
        disp(['- ' num2str(rangeK(k)) ' FC states'])
        for c = 1:rangeK(k)
            cond_pair = 1;
            for cond1 = 1:n_Cond-1
                for cond2 = cond1+1:n_Cond
                    a = squeeze(P(Index_Conditions == cond1,k,c));  % Vector containing Prob of c in condition 1
                    b = squeeze(P(Index_Conditions == cond2,k,c));  % Vector containing Prob of c in condition 2

                    % data points of prob of occurence for each clustering
                    % solution from different conditions
                    data_vec = cat(1,a,b);
                    % column vector with index of groups
                    groups = cat(1,zeros(numel(a),1),ones(numel(b),1));
                    % perform Levene's Test of equality of variances (store p-value)
                    levene_pval(cond_pair,k,c) = vartestn(data_vec,groups,'TestType','LeveneAbsolute','Display','off');
                    cond_pair = cond_pair + 1;
                end
            end
        end
    end
end

%% PERMUTATION STATISTICS WITH WITHIN BOOTSTRAP SAMPLES

disp('The hypothesis tests take a considerable amount of time to run.')
disp('Testing differences in fractional occupancy:')

% Store p-values and effect size for K-means clustering solutions
P_pval = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));
P_pval2sided = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));
effectsize = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));

for k = 1:length(rangeK)
    disp(['- ' num2str(rangeK(k)) ' FC states'])
    for c = 1:rangeK(k)
        cond_pair = 1;
        for cond1 = 1:n_Cond-1
            for cond2 = cond1+1:n_Cond
                a = squeeze(P(Index_Conditions == cond1,k,c))';  % Vector containing Prob of c in condition 1
                b = squeeze(P(Index_Conditions == cond2,k,c))';  % Vector containing Prob of c in condition 2
                
                if pair == 1
                    % disp([cond{cond1} ' vs ' cond{cond2} ' (using paired-sample t-test statistic)'])
                    stats = bootstrap_within_permutation_paired_samples([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],...,
                                                                        n_permutations,n_bootstraps,0.05);
                else
                    if levene_pval(cond_pair,k,c) < 0.05 % null hypothesis of homogeneity of variances rejected
                        % disp([cond{cond1} ' vs ' cond{cond2} ' (using Welch''s t-test statistic)'])
                        stats = bootstrap_within_permutation_ttest2([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],n_permutations,...,
                                                                     n_bootstraps,0.05,'welchtest');
                    else % null hypothesis of homegeneity of variances not rejected
                        % disp([cond{cond1} ' vs ' cond{cond2} ' (using t-test statistic)'])
                        stats = bootstrap_within_permutation_ttest2([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],n_permutations,...,
                                                                     n_bootstraps,0.05,'ttest');
                    end
                end

                % p-values and effect size
                P_pval(cond_pair,k,c) = min(stats.pvals);
                P_pval2sided(cond_pair,k,c) = 2*min(stats.pvals);
                effectsize(cond_pair,k,c) = stats.eff;
                cond_pair = cond_pair + 1;
            end
        end
    end
end

% Name of the file to save output
save_file = 'LEiDA_Stats_FracOccup.mat';

% Save K-means clustering solutions results:
save([data_dir '/' save_file],'P','P_pval','P_pval2sided', 'effectsize', 'levene_pval',...,
                              'cond','rangeK','file_cluster','file_V1','Index_Conditions','pair')

disp(['Fractional occupancy values and hypothesis tests results saved successfully as ' save_file])