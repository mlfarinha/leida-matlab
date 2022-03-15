function LEiDA_stats_TransitionMatrix(data_dir,bestK,pair)
%
% For the optimalK compute the transition matrix for each participant.
% Perform hypothesis tests to check for differences in the state-to-state
% transition probabilities between conditions.
%
% INPUT:
% data_dir       directory where the LEiDA results are saved
% bestK          optimal K defined by the user
% pair           0, different subjects in each condition (default); 1, same
%                subjects across conditions
%
% OUTPUT:
% TM             transition matrix for optimal K for all participants
% TMnorm         transition probability matrix (TPM) for optimal K for all
%                participants
% TM_pval2sided  two-sided p-value from testing whether the mean transition
%                probability of a given FC state differs between conditions      
% effectsize    Hedge's effect size
% levene_pval   p-value obtained from computing the Levene's test on the
%               state-to-state transition probability values
%
% Author: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% Default number of permutations
n_permutations = 500;
% Default number of bootstrap samples within each permutation sample
n_bootstraps = 10;

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';

% Load required data:
load([data_dir file_V1], 'Time_sessions');
load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
load([data_dir file_P], 'cond', 'P', 'Index_Conditions');

N_scans = max(Time_sessions);

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Code below computes the fractional occupancy and runs hypothesis tests
disp('%%%%%%%%%%%%%%%%%%%%%%% TRANSITION PROBABILITIES %%%%%%%%%%%%%%%%%%%%%%%')

%% CALCULATE THE TRANSITION MATRIX FOR THE OPTIMAL K

% Matrix with probability of occurrence of each FC state for chosen K
P_K = squeeze(P(:,bestK-1,1:bestK));

% Transition probability matrix
TM = zeros(N_scans, bestK, bestK);
TMnorm = zeros(N_scans, bestK, bestK);

for s = 1:N_scans

    T = Time_sessions == s;
    Ctime = Kmeans_results{rangeK==bestK}.IDX(T);
    % alpha corresponds to departure state;
    for alpha = 1:bestK
        % beta corresponds to arrival state;
        for beta = 1:bestK
            % count number transition from alpha to beta
            alpha2beta = 0;
            for i = 1:length(Ctime)-1
                % if alpha = beta then we have 1 transition
                if isequal(Ctime(i),alpha) && isequal(Ctime(i+1),beta)
                    alpha2beta = alpha2beta + 1;
                end
            end
            TM(s,alpha,beta) = alpha2beta;
        end
    end
    % normalise by T-1 as in Vohryzek et al., (2020)
    TM(s,:,:) = TM(s,:,:)./(size(Ctime,2)-1);
    TMnorm(s,:,:) = squeeze(TM(s,:,:))./P_K(s,:)';
end

% Set all Nan values to zero
% If a FC state has probability of occurence 0 then there should be a
% probability of 0 to transition to that state
% numNaN = sum(isnan(TMnorm),'all');
% TMnorm(isnan(TMnorm)) = 0;

clear Kmeans_results

%% CHECK HOMOGENEITY OF VARIANCES USING LEVENE'S TEST

if pair == 0
    
    disp('Running Levene''s test of equality of variances:')
    % matrix to store p-value of Levene's test of equality of variances
    levene_pval = zeros(n_Cond*(n_Cond-1)/2,bestK,bestK);

    for c_out = 1:bestK
        for c_in = 1:bestK
            disp(['- Transition from FC state ' num2str(c_out) ' to FC state ' num2str(c_in)])
            cond_pair = 1;
            for cond1 = 1:n_Cond-1
                for cond2 = cond1+1:n_Cond
                    a = squeeze(TMnorm(Index_Conditions == cond1,c_out,c_in));
                    b = squeeze(TMnorm(Index_Conditions == cond2,c_out,c_in));
                    
                    data_vec = cat(1,a,b);
                    % column vector with index of groups, SZ = 0 and HC = 1
                    groups = cat(1,zeros(numel(a),1),ones(numel(b),1));
                    % perform Levene Test of equality of variances (store p-value)
                    levene_pval(cond_pair,c_out,c_in) = vartestn(data_vec,groups,'TestType','LeveneAbsolute','Display','off');
                    cond_pair = cond_pair + 1;
                end
            end
        end
    end
end

%%  PERMUTATION STATISTICS WITH WITHIN BOOTSTRAP SAMPLES

disp('The hypothesis tests take a considerable amount of time to run')
disp('Testing differences in transition probabilities:')

% Store p-values and effect size
TM_pval = zeros(n_Cond*(n_Cond-1)/2,bestK,bestK);
TM_pval2sided = zeros(n_Cond*(n_Cond-1)/2,bestK,bestK);
effectsize = zeros(n_Cond*(n_Cond-1)/2,bestK,bestK);

for c_out = 1:bestK 
    for c_in = 1:bestK
        disp(['- Transition from FC state ' num2str(c_out) ' to FC state ' num2str(c_in)])
        cond_pair = 1;
        for cond1 = 1:n_Cond-1
            for cond2 = cond1+1:n_Cond
                a = squeeze(TMnorm(Index_Conditions == cond1,c_out,c_in))';
                a = rmmissing(a);
                b = squeeze(TMnorm(Index_Conditions == cond2,c_out,c_in))';  % Vector containing Prob of c in condition 2
                b = rmmissing(b);

                if pair == 1
                    % disp([cond{cond1} ' vs ' cond{cond2} ' (using paired-sample t-test statistic)'])
                    stats = bootstrap_within_permutation_paired_samples([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],...,
                                                                        n_permutations,n_bootstraps,0.05);
                else      
                    if levene_pval(cond_pair,c_out,c_in) < 0.05 % null hypothesis of homogeneity of variances rejected
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
                TM_pval(cond_pair,c_out,c_in) = min(stats.pvals);
                TM_pval2sided(cond_pair,c_out,c_in) = 2*min(stats.pvals);
                effectsize(cond_pair,c_out,c_in) = stats.eff;
                cond_pair = cond_pair + 1;
            end
        end
    end
end

% Name of the file to save output
save_file = ['LEiDA_Stats_TransitionMatrix_K' num2str(bestK) '.mat'];

save([data_dir '/' save_file],'TM','TMnorm','TM_pval','TM_pval2sided','effectsize','levene_pval',...,
                              'cond','rangeK','file_V1','file_cluster','file_P','Index_Conditions')

disp(['Transition probability values and hypothesis tests results saved successfully as ' save_file])