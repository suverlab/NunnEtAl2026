%% 11/25/25 OMN
% creates APN2 Vm |AUC| plot
    % combines LVm and RVm values for each condition (MAT FILE)
        % control no glue, no glue, control glue, glue
    % the original sign (positive VS negative) of each dot is represented by its respective color 
%% to call:
% plot_auc_abs_combined_signColored('control_noGlueQ_Vm_AUCs.mat', 'noGlueQ_Vm_AUCs.mat', 'control_glueQ_Vm_AUCs.mat', 'glueQ_Vm_AUCs.mat')

%%
function plot_auc_abs_combined_signColored(file1, file2, file3, file4)

    files = {file1, file2, file3, file4};
    labels = {'Ctrl NoGlueQ', 'NoGlueQ', 'Ctrl GlueQ', 'GlueQ'};
    groupColors = [0 0 1; 0 0.5 1; 1 0 0; 1 0.5 0]; % used for mean+SEM only
    jitter_amount = 0.1;

    % Colors for the dots based on SIGN
    posColor = [0    0.4  1];    % positive original AUC → blue
    negColor = [1    0.2  0.2];  % negative original AUC → red

    figure; hold on;

    for i = 1:4
        data = load(files{i});

        % -------- Combine L & R AUCs into a single vector --------
        auc_raw = [ ...
            data.Lvm_pos_auc, data.Rvm_pos_auc, ...
            data.Lvm_neg_auc, data.Rvm_neg_auc ...
        ];

        % -------- Absolute values for plotting height --------
        auc_abs = abs(auc_raw);

        % -------- Jitter x-positions --------
        x_jitter = i + jitter_amount*randn(size(auc_raw));

        % -------- Scatter each point with sign-dependent color --------
        for k = 1:length(auc_raw)
            if auc_raw(k) >= 0
                dotColor = posColor;
            else
                dotColor = negColor;
            end

            scatter(x_jitter(k), auc_abs(k), 50, ...
                'MarkerFaceColor', dotColor, ...
                'MarkerFaceAlpha', 0.7, ...
                'MarkerEdgeColor', 'none');
        end

        % -------- Mean + SEM of ABSOLUTE values --------
        m = mean(auc_abs);
        sem = std(auc_abs) / sqrt(numel(auc_abs));

        % Horizontal mean line
        plot([i-0.2 i+0.2], [m m], 'k', 'LineWidth', 2);

        % Vertical SEM line + caps
        line([i i], [m-sem m+sem], 'Color','k','LineWidth',2);
        line([i-0.05 i+0.05], [m-sem m-sem], 'Color','k','LineWidth',2);
        line([i-0.05 i+0.05], [m+sem m+sem], 'Color','k','LineWidth',2);

        % n count
        text(i, m + sem + 0.02, sprintf('n=%d', numel(auc_abs)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontWeight','bold', 'FontSize', 9);
    end

    % -------- Axes + Labels --------
    xlim([0.5 4.5]);
    xticks(1:4);
    xticklabels(labels);
    ylabel('|AUC| (absolute)');
    title('Absolute LVm + RVm AUCs (Dot Color Indicates Original Sign)');

    % -------- Legend --------
        % Create invisible marker handles with the correct colors and marker
    hPos = plot(NaN, NaN, 'o', 'MarkerFaceColor', posColor, ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    hNeg = plot(NaN, NaN, 'o', 'MarkerFaceColor', negColor, ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 8);

    legend([hPos, hNeg], {'Originally positive','Originally negative'}, ...
           'Location','northwest', 'Box','off');

end
