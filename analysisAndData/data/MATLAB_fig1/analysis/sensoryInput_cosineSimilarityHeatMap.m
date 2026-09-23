%% APN2 Input Cosine Similarity with Weights (Using Synapse Counts)

apn2IDs_sorted = [apnMirrorPairs(:,1);apnMirrorPairs(:,2)];
apn2IDs_sorted(strcmp(apn2IDs_sorted,"no match")) = [];

% Create weighted matrix
weighted_matrix = zeros(numel(apn2IDs_sorted), numel(sInputID));

% Populate the weighted matrix using synapse counts
for ii = 1:numel(apn2IDs_sorted)
    currentAPN2 = apn2IDs_sorted(ii);
    APN2conn = apnInputConns(apnInputConns.post_root_id == currentAPN2, :);

    for jj = 1:height(APN2conn)
        neuron_idx = find(strcmp(sInputID, APN2conn.pre_root_id{jj}));
        weighted_matrix(ii, neuron_idx) = weighted_matrix(ii, neuron_idx) + APN2conn.syn_count(jj);
    end
end

% Compute cosine similarity
cosine_similarity_weighted = zeros(numel(apn2IDs_sorted));
for ii = 1:numel(apn2IDs_sorted)
    for jj = 1:numel(apn2IDs_sorted)
        if norm(weighted_matrix(ii, :)) > 0 && norm(weighted_matrix(jj, :)) > 0
            cosine_similarity_weighted(ii, jj) = dot(weighted_matrix(ii, :), weighted_matrix(jj, :)) / ...
                                               (norm(weighted_matrix(ii, :)) * norm(weighted_matrix(jj, :)));
        else
            cosine_similarity_weighted(ii, jj) = 0; % Handle zero norm cases
        end
    end
end

% Visualize cosine similarity matrix
figCosine = figure('Color', 'white', 'Position', [500, 50, 700, 650]);
imagesc(cosine_similarity_weighted);
colormap(flipud(gray));  % Apply grayscale colormap
colorbar;
xticks(1:numel(apn2IDs_sorted));
yticks(1:numel(apn2IDs_sorted));
xticklabels(apn2IDs_sorted);
yticklabels(apn2IDs_sorted);
title('Weighted Cosine Similarity of APN2s and Sensory Inputs');
set(gca, 'XTickLabelRotation', 90);

% Add cosine similarity values to the cells
for ii = 1:numel(apn2IDs_sorted)
    for j = 1:numel(apn2IDs_sorted)
        value = cosine_similarity_weighted(ii, j);
        text_color = 'w';  % default black
        if value < 0.5     % if background is dark, switch to white
            text_color = 'k';
        end
        text(j, ii, sprintf('%.2f', value), ...
            'HorizontalAlignment', 'center', 'Color', text_color, 'FontSize', 10);
    end
end

% save figures
cd(savedFiguresPath)

saveas(figCosine,todaysDate+"_apnCosineSimilarityMatrix",'fig')
saveas(figCosine,todaysDate+"_apnCosineSimilarityMatrix",'png')
saveas(figCosine,todaysDate+"_apnCosineSimilarityMatrix",'pdf')

% Perform Pairwise Comparisons using Mann-Whitney U Test
pairwise_p_values = zeros(numel(apn2IDs_sorted));
for ii = 1:numel(apn2IDs_sorted)
    for jj = ii+1:numel(apn2IDs_sorted)
        p_mannwhitney = ranksum(cosine_similarity_weighted(ii, :), cosine_similarity_weighted(j, :));
        pairwise_p_values(ii, jj) = p_mannwhitney;
        pairwise_p_values(jj, ii) = p_mannwhitney; % Symmetric matrix
    end
end

% Display pairwise p-values
disp('Pairwise Mann-Whitney U Test p-values:');
disp(array2table(pairwise_p_values, 'VariableNames',apn2IDs_sorted, 'RowNames', apn2IDs_sorted));
