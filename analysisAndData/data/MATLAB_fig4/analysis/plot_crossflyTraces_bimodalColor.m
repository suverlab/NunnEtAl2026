%% OMN 11/24/25
    % plots the average Laa, Raa, Lvm, Rvm across all exps
        %  !!! with negative Vm responders plotted in one color and positive in another !!!
        % w/global average plotted on top
%% instructions:
% 1) data = load('yourMAT.mat') % make sure you youve labeled ROIs with GUI
% 2) groupedTrials = group_trials_by_date_exp(data)
% 3) comboAvg = compute_combo_averages(data, groupedTrials, 20)
function plot_crossflyAvgs_bimodalColor(comboAvg)
    
    figure;

    % --- Create 2x2 subplots ---
    ax = gobjects(2,2);
    ax(1,1) = subplot(2,2,1); hold on; title('L Antenna'); xlabel('samples'); ylabel('Δant position');
    ax(1,2) = subplot(2,2,2); hold on; title('R Antenna'); xlabel('samples'); ylabel('Δant position');
    ax(2,1) = subplot(2,2,3); hold on; title('L Vm (pos+neg)'); xlabel('samples'); ylabel('Vm (mV)');
    ax(2,2) = subplot(2,2,4); hold on; title('R Vm (pos+neg)'); xlabel('samples'); ylabel('Vm (mV)');


    %% color scheme
    posColor = '#56A5FF';  % light blue
    negColor = '#1A4DEC';  % dark blue

    % Counts: row1=antenna, row2=Vm positive, row3=Vm negative
    counts = zeros(3,2);

    % Storage for averages
    Laa_all = {}; Raa_all = {};
    Lvm_pos_all = {}; Lvm_neg_all = {};
    Rvm_pos_all = {}; Rvm_neg_all = {};

    % AUC storage
    Lvm_pos_auc = [];
    Lvm_neg_auc = [];
    Rvm_pos_auc = [];
    Rvm_neg_auc = [];


    % ==============================
    %        LOOP EXPERIMENTS
    % ==============================
    for e = 1:numel(comboAvg)

        % Determine how many trials to loop through (max of antenna/vm)
        L_n = max([numel(comboAvg(e).Laa), numel(comboAvg(e).Lvm)]);
        R_n = max([numel(comboAvg(e).Raa), numel(comboAvg(e).Rvm)]);

        % ------------------------------
        %            LEFT SIDE
        % ------------------------------
        for t = 1:L_n

            % --- Determine trace color from Vm trial ---
            if t <= numel(comboAvg(e).Lvm)
                vmTrace = comboAvg(e).Lvm{t};
                auc = trapz(vmTrace);

                if auc >= 0
                    colorToUse = posColor;
                    % Store positive Vm
                    Lvm_pos_all{end+1} = vmTrace;
                    Lvm_pos_auc(end+1) = auc;
                    counts(2,1) = counts(2,1) + 1;
                else
                    colorToUse = negColor;
                    % Store negative Vm
                    Lvm_neg_all{end+1} = vmTrace;
                    Lvm_neg_auc(end+1) = auc;
                    counts(3,1) = counts(3,1) + 1;
                end

                % Plot Vm trace
                plot(ax(2,1), vmTrace, 'Color', colorToUse, 'LineWidth', 1.5);
            else
                % If no Vm trace exists, assume positive color
                colorToUse = posColor;
            end

            % --- Plot antenna trial with same color (if exists) ---
            if t <= numel(comboAvg(e).Laa)
                plot(ax(1,1), comboAvg(e).Laa{t}, 'Color', colorToUse, 'LineWidth', 1.5);
                Laa_all{end+1} = comboAvg(e).Laa{t};
                counts(1,1) = counts(1,1) + 1;
            end
        end


        % ------------------------------
        %            RIGHT SIDE
        % ------------------------------
        for t = 1:R_n

            % Determine Vm AUC → choose color
            if t <= numel(comboAvg(e).Rvm)
                vmTrace = comboAvg(e).Rvm{t};
                auc = trapz(vmTrace);

                if auc >= 0
                    colorToUse = posColor;
                    Rvm_pos_all{end+1} = vmTrace;
                    Rvm_pos_auc(end+1) = auc;
                    counts(2,2) = counts(2,2) + 1;
                else
                    colorToUse = negColor;
                    Rvm_neg_all{end+1} = vmTrace;
                    Rvm_neg_auc(end+1) = auc;
                    counts(3,2) = counts(3,2) + 1;
                end

                % Plot Vm trace
                plot(ax(2,2), vmTrace, 'Color', colorToUse, 'LineWidth', 1.5);
            else
                colorToUse = posColor; % fallback
            end

            % Plot antenna corresponding color
            if t <= numel(comboAvg(e).Raa)
                plot(ax(1,2), comboAvg(e).Raa{t}, 'Color', colorToUse, 'LineWidth', 1.5);
                Raa_all{end+1} = comboAvg(e).Raa{t};
                counts(1,2) = counts(1,2) + 1;
            end
        end
    end


    % ==============================
    %        TRIM AUCs TO 150
    % ==============================
    maxPoints = 150;
    Lvm_pos_auc = Lvm_pos_auc(1:min(end,maxPoints));
    Lvm_neg_auc = Lvm_neg_auc(1:min(end,maxPoints));
    Rvm_pos_auc = Rvm_pos_auc(1:min(end,maxPoints));
    Rvm_neg_auc = Rvm_neg_auc(1:min(end,maxPoints));

    save('Vm_AUCs.mat', 'Lvm_pos_auc','Lvm_neg_auc','Rvm_pos_auc','Rvm_neg_auc');


    % ==============================
    %        COUNTS TEXT
    % ==============================

    % Top-left: Left antenna count (counts row 1, col 1)
    text(ax(1,1), 0.02, 0.95, sprintf('n = %d', counts(1,1)), ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold');
    
    % Top-right: Right antenna count (counts row 1, col 2)
    text(ax(1,2), 0.98, 0.95, sprintf('n = %d', counts(1,2)), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold');
    
    % Bottom-left (L Vm): show Positive and Negative counts separately
    text(ax(2,1), 0.02, 0.95, sprintf('pos n = %d', counts(2,1)), ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold');
    text(ax(2,1), 0.98, 0.95, sprintf('neg n = %d', counts(3,1)), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold');
    
    % Bottom-right (R Vm): show Positive and Negative counts separately
    text(ax(2,2), 0.02, 0.95, sprintf('pos n = %d', counts(2,2)), ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold');
    text(ax(2,2), 0.98, 0.95, sprintf('neg n = %d', counts(3,2)), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold');

    % ==============================
    %        AXIS + TICKS
    % ==============================
    % Antenna
    ylim(ax(1,1), [-12 35]);
    ylim(ax(1,2), [-12 35]);
    yticks(ax(1,1), [0 10 20 30]);
    yticks(ax(1,2), [0 10 20 30]);

    % Vm combined
    ylim(ax(2,1), [-12 7]);
    ylim(ax(2,2), [-12 7]);
    yticks(ax(2,1), [-9 -6 -3 0 3 6]);
    yticks(ax(2,2), [-9 -6 -3 0 3 6]);

    for r = 1:2
        for c = 1:2
            xlim(ax(r,c), [0 150]);
            xticks(ax(r,c), [0 50 100 150]);
        end
    end

    sgtitle('All Trials — Antenna and Vm with Counts (Combined Vm)');

end

