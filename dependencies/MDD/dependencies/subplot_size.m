function [subplot_rows,subplot_cols]=subplot_size(no_plots)

subplot_cols=ceil(sqrt(no_plots));
subplot_rows=ceil(no_plots/subplot_cols);