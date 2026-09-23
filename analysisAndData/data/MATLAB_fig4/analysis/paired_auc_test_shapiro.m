%% 11/24/25 OMN
% performs normality tests on all data groups
    % assigns respective stats test for paired and independent comparisons
%% to call: 
% paired_auc_test_shapiro('control_noGlueQ_Vm_AUCs.mat', 'noGlueQ_Vm_AUCs.mat', 'control_glueQ_Vm_AUCs.mat', 'glueQ_Vm_AUCs.mat')

function paired_auc_test_shapiro(ctrlNoGlueFile, noGlueFile, ctrlGlueFile, glueFile)
    % Load MAT files
    ctrlNoGlue = load(ctrlNoGlueFile);
    noGlue = load(noGlueFile);
    ctrlGlue = load(ctrlGlueFile);
    glue = load(glueFile);

    % Combine L+R AUCs
    auc_ctrlNoGlue = [ctrlNoGlue.Lvm_pos_auc, ctrlNoGlue.Rvm_pos_auc, ...
                      ctrlNoGlue.Lvm_neg_auc, ctrlNoGlue.Rvm_neg_auc];
    auc_noGlue = [noGlue.Lvm_pos_auc, noGlue.Rvm_pos_auc, ...
                  noGlue.Lvm_neg_auc, noGlue.Rvm_neg_auc];
    auc_ctrlGlue = [ctrlGlue.Lvm_pos_auc, ctrlGlue.Rvm_pos_auc, ...
                    ctrlGlue.Lvm_neg_auc, ctrlGlue.Rvm_neg_auc];
    auc_glue = [glue.Lvm_pos_auc, glue.Rvm_pos_auc, ...
                glue.Lvm_neg_auc, glue.Rvm_neg_auc];

    % Trim to same length for paired test
    n1 = min(length(auc_ctrlNoGlue), length(auc_noGlue));
    n2 = min(length(auc_ctrlGlue), length(auc_glue));
    auc_ctrlNoGlue = auc_ctrlNoGlue(1:n1);
    auc_noGlue = auc_noGlue(1:n1);
    auc_ctrlGlue = auc_ctrlGlue(1:n2);
    auc_glue = auc_glue(1:n2);

    %% if you wants stats on ABS value!!
    auc_ctrlNoGlue = abs(auc_ctrlNoGlue);
    auc_noGlue = abs(auc_noGlue);
    auc_ctrlGlue = abs(auc_ctrlGlue);
    auc_glue = abs(auc_glue);

    % --- Test normality using Shapiro-Wilk ---
    alpha = 0.05;
    [h1_ctrl, p1_ctrl, W1_ctrl] = swtest(auc_ctrlNoGlue, alpha);  % ctrl no-glue
    [h1_test,  p1_test,  W1_test]  = swtest(auc_noGlue, alpha);

    [h2_ctrl, p2_ctrl, W2_ctrl] = swtest(auc_ctrlGlue, alpha);
    [h2_test,  p2_test,  W2_test]  = swtest(auc_glue, alpha);

    fprintf('Shapiro-Wilk test (NoGlue): Ctrl: p=%.3f, W=%.3f, h=%d; Test: p=%.3f, W=%.3f, h=%d\n', ...
        p1_ctrl, W1_ctrl, h1_ctrl, p1_test, W1_test, h1_test);
    fprintf('Shapiro-Wilk test (Glue): Ctrl: p=%.3f, W=%.3f, h=%d; Test: p=%.3f, W=%.3f, h=%d\n', ...
        p2_ctrl, W2_ctrl, h2_ctrl, p2_test, W2_test, h2_test);

    % --- Paired comparison ---
    if (h1_ctrl == 0) && (h1_test == 0)
        [~, p1, ~, stats1] = ttest(auc_ctrlNoGlue, auc_noGlue);
        fprintf('Paired t-test (NoGlue): t(%d)=%.2f, p=%.4f\n', stats1.df, stats1.tstat, p1);
    else
        [p1, ~, stats1] = signrank(auc_ctrlNoGlue, auc_noGlue);
        fprintf('Wilcoxon signed-rank (NoGlue): p=%.4f\n', p1);
    end

    if (h2_ctrl == 0) && (h2_test == 0)
        [~, p2, ~, stats2] = ttest(auc_ctrlGlue, auc_glue);
        fprintf('Paired t-test (Glue): t(%d)=%.2f, p=%.4f\n', stats2.df, stats2.tstat, p2);
    else
        [p2, ~, stats2] = signrank(auc_ctrlGlue, auc_glue);
        fprintf('Wilcoxon signed-rank (Glue): p=%.4f\n', p2);
    end

    % --- Independent Tests (NoGlue vs Glue) ---
    fprintf('\n===== Independent comparisons (Ctrl NoGlue vs Ctrl Glue) =====\n');
    % Normality
    h_ctrlNG = swtest(auc_ctrlNoGlue, alpha);
    h_ctrlG  = swtest(auc_ctrlGlue, alpha);

    if h_ctrlNG == 0 && h_ctrlG == 0
        [~, p_ind1, ~, stats_ind1] = ttest2(auc_ctrlNoGlue, auc_ctrlGlue);
        fprintf('Independent t-test: t(%d)=%.2f, p=%.4f\n', stats_ind1.df, stats_ind1.tstat, p_ind1);
    else
        p_ind1 = ranksum(auc_ctrlNoGlue, auc_ctrlGlue);
        fprintf('Mann-Whitney U (ranksum): p=%.4f\n', p_ind1);
    end

    fprintf('\n===== Independent comparisons (Test NoGlue vs Test Glue) =====\n');
    h_testNG = swtest(auc_noGlue, alpha);
    h_testG  = swtest(auc_glue, alpha);

    if h_testNG == 0 && h_testG == 0
        [~, p_ind2, ~, stats_ind2] = ttest2(auc_noGlue, auc_glue);
        fprintf('Independent t-test: t(%d)=%.2f, p=%.4f\n', stats_ind2.df, stats_ind2.tstat, p_ind2);
    else
        p_ind2 = ranksum(auc_noGlue, auc_glue);
        fprintf('Mann-Whitney U (ranksum): p=%.4f\n', p_ind2);
    end

end
