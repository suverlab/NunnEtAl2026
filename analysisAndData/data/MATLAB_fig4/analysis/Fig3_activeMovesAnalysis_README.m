%% README
% Nunn et al., 2026


%% Figure 3 A-B: Active movement plotting
%% 
% 1. download MAT files containing Vm and antennal data correlated to active movement indicies:
% '200ms_LandR_activeMoves_glue_Qtrials.mat'
% '200ms_LandR_activeMoves_noGlue_Qtrials.mat'
% 'control_LandR_activeMoves_glue_Qtrials.mat'
% 'control_LandR_activeMoves_noGlue_Qtrials.mat'

%%
% 2. in MATLAB load data file containing active movement indicies
% ex) data = load('200ms_LandR_activeMoves_glue_Qtrials.mat')

%%
% 3. combine experiment date and experiment number to have
% unique identifiers for each cell (N)
%  ex) groupedTrials = group_trials_by_date_exp(data))

%%
% 4. calculate the average Laa, Raa, Lvm, Rvm across trials for each experiment
% ex) comboAvg = compute_combo_averages(data,groupedTrials,buffer)
    % buffer = 20 (amount of samples to plot before and after ROI)

%%
% 5. plots Laa, Raa, Lvm, Rvm traces, with the group average overlaid
% B) colors of aa match vm, all hyperpolarizing trials are one color and depolarizing trials another color
    % ex) plot_crossflyTraces_bimodalColor(comboAvg)


%% Figure 3C: Active movement quantification
% 6. Download MAT files containing calculated area under the curve for each condition:
% control_noGlueQ_Vm_AUCs.mat
% control_glueQ_Vm_AUCs.mat
% noGlueQ_Vm_AUCs.mat
% glueQ_Vm_AUCs.mat

%%
% 7. plot quantification (APN2 Vm |AUC|)
% combines LVm (ipsi) and RVm (contra) values for each condition: (control no glue, no glue, control glue, glue)
% the original sign (positive/depolarizing VS negative/hyperpolarizing) of each dot is represented by its respective color 
% ex) plot_auc_abs_combined_signColored('control_noGlueQ_Vm_AUCs.mat', 'noGlueQ_Vm_AUCs.mat', 'control_glueQ_Vm_AUCs.mat', 'glueQ_Vm_AUCs.mat')

%%
% 8. compute stats
% ex) paired_auc_test_shapiro('control_noGlueQ_Vm_AUCs.mat', 'noGlueQ_Vm_AUCs.mat', 'control_glueQ_Vm_AUCs.mat', 'glueQ_Vm_AUCs.mat')
     % need to downlaod the swtest package from MATLAB community center,prior


