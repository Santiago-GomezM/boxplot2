# boxplot2

boxplot2 plots a 2D boxplot for the given data and categories (species).

The function expects the input 'data' to be a matrix with dimensions m-by-2 or 2-by-n,
where m or n can be any number but one of the dimensions must be equal to 2.
The 'species' input must be a categorical array that is a column vector.
The function will plot the boxplots for the two variables (columns) in 'data', color-coded by species, and will mark the outliers.
