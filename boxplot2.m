function boxplot2(data, species)

    % boxplot2 plots a 2D boxplot for the given data and categories (species).
    % The function expects the input 'data' to be a matrix with dimensions m-by-2 or 2-by-n,
    % where m or n can be any number but one of the dimensions must be equal to 2.
    % The 'species' input must be a categorical array that is a column vector.
    % The function will plot the boxplots for the two variables (columns) in 'data', 
    % color-coded by species, and will mark the outliers.

    % Author: Santiago Gómez
    % Last Update: 07/11/2024
    % Created: 07/11/2024
    
    % Check if 'data' has one dimension equal to 2
    if size(data, 2) == 2
        % Data is already in the correct form (m-by-2)
    elseif size(data, 1) == 2
        % Transpose data to make sure the second dimension is 2 (2-by-n)
        data = data';
    else
        error('Data must have one dimension equal to 2.');
    end
    
    % Check if 'species' is a categorical array and reshape it to a column if necessary
    if ~iscategorical(species)
        error('The species array must be of type categorical.');
    end
    
    if size(species, 2) ~= 1
        species = species(:);  % Reshape species to be a column vector
        %warning('The species array has been reshaped to a column.');
    end
    
    % Calculate boxplot statistics manually
    stats = struct();
    speciesNames = unique(species);  % Get unique species
    
    for i = 1:numel(speciesNames)
        idx = species == speciesNames(i);  % Find indices for the current species
        stats(i).category = speciesNames(i);
        
        % Calculate statistics for Sepal Length (x)
        stats(i).x_q1 = quantile(data(idx, 1), 0.25);
        stats(i).x_median = median(data(idx, 1));
        stats(i).x_mean = mean(data(idx, 1));
        stats(i).x_q3 = quantile(data(idx, 1), 0.75);
        
        % Calculate statistics for Sepal Width (y)
        stats(i).y_q1 = quantile(data(idx, 2), 0.25);
        stats(i).y_median = median(data(idx, 2));
        stats(i).y_mean = mean(data(idx, 2));
        stats(i).y_q3 = quantile(data(idx, 2), 0.75);
        
        % Calculate outliers manually using IQR method
        x_iqr = stats(i).x_q3 - stats(i).x_q1;
        y_iqr = stats(i).y_q3 - stats(i).y_q1;
        stats(i).y_min = stats(i).y_q1 - 1.5 * y_iqr;
        stats(i).y_max = stats(i).y_q3 + 1.5 * y_iqr;
        stats(i).x_min = stats(i).x_q1 - 1.5 * x_iqr;
        stats(i).x_max = stats(i).x_q3 + 1.5 * x_iqr;

        stats(i).x_outliers = data(idx, 1) < (stats(i).x_q1 - 1.5 * x_iqr) | ...
                              data(idx, 1) > (stats(i).x_q3 + 1.5 * x_iqr);
        stats(i).y_outliers = data(idx, 2) < (stats(i).y_q1 - 1.5 * y_iqr) | ...
                              data(idx, 2) > (stats(i).y_q3 + 1.5 * y_iqr);
    end

    % Plot the boxplots with outliers
    figure;
    hold on;
    colors = lines(numel(speciesNames));  % Generate colors for each species
    legendEntries = {};  % List to hold legend entries
    patchHandles = zeros(numel(speciesNames), 1);  % List to hold patch handles

    for i = 1:numel(stats)
        % Draw the box (rectangle) using patch
        x_coords = [stats(i).x_q1, stats(i).x_q1, stats(i).x_q3, stats(i).x_q3];
        y_coords = [stats(i).y_q1, stats(i).y_q3, stats(i).y_q3, stats(i).y_q1];
        
        % Use patch to draw the rectangle with transparency
        patchHandles(i) = patch(x_coords, y_coords, colors(i,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 0.75);
        
        % Add the species category to the legend
        legendEntries{i} = char(stats(i).category);
        
        % Plot the whiskers for x dimension
        plot([stats(i).x_min stats(i).x_max], [stats(i).y_median stats(i).y_median], 'k-', 'LineWidth', 0.75, 'Color', colors(i,:) * 0.7);
        plot([stats(i).x_min stats(i).x_min], [stats(i).y_q1 stats(i).y_q3], 'k-', 'LineWidth', 0.75, 'Color', colors(i,:) * 0.7);
        plot([stats(i).x_max stats(i).x_max], [stats(i).y_q1 stats(i).y_q3], 'k-', 'LineWidth', 0.75, 'Color', colors(i,:) * 0.7);
        
        % Plot the whiskers for y dimension
        plot([stats(i).x_median stats(i).x_median], [stats(i).y_min stats(i).y_max], 'k-', 'LineWidth', 0.75, 'Color', colors(i,:) * 0.7);
        plot([stats(i).x_q1 stats(i).x_q3], [stats(i).y_min stats(i).y_min], 'k-', 'LineWidth', 0.75, 'Color', colors(i,:) * 0.7);
        plot([stats(i).x_q1 stats(i).x_q3], [stats(i).y_max stats(i).y_max], 'k-', 'LineWidth', 0.75, 'Color', colors(i,:) * 0.7);

        % Plot the outliers as crosses (x)
        x_outliers = data(stats(i).x_outliers,1);
        y_outliers = data(stats(i).y_outliers,2);

        x_outliers = [x_outliers ones(numel(x_outliers),1)*stats(i).y_median];
        y_outliers = [ones(numel(y_outliers),1)*stats(i).x_median y_outliers];
                
        % Plot the outliers with crosses
        scatter(x_outliers(:,1), x_outliers(:,2), 50, 'Marker', 'x', 'MarkerEdgeColor', colors(i,:)*0.7, 'LineWidth', 0.75);
        scatter(y_outliers(:,1), y_outliers(:,2), 50, 'Marker', 'x', 'MarkerEdgeColor', colors(i,:)*0.7, 'LineWidth', 0.75);

        % Plot the mean with circle
        scatter(stats(i).x_mean, stats(i).y_mean, 50, 'Marker', 'o', 'MarkerEdgeColor', colors(i,:)*0.7, 'LineWidth', 0.75);
    end

    % Add the legend with species categories and colors
    legend(patchHandles, legendEntries, 'Location', 'eastoutside');

    % Labels and title
    hold off;
    grid on;
    box on;

end

