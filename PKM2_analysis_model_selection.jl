using CSV, DataFrames, CairoMakie, Statistics

# res_df = CSV.read("Cluster_results/100823_PKM2_forward_subset_selection_train_test_loss_for_all_figs_removed_params_top_10_percent.csv", DataFrame)
# res_df = CSV.read("Cluster_results/100823_PKM2_reverse_subset_selection_train_test_loss_for_all_figs_removed_params_top_10_percent.csv", DataFrame)
# res_df.num_params =
#     [23 - sum(values(eval(Meta.parse(param_subset))) .> 0) for param_subset in res_df.param_subset]
# res_df = CSV.read("Cluster_results/100923_PKM2_forward_subset_selection_train_test_loss_for_all_figs_5_repeats.csv", DataFrame)
# res_df = CSV.read("Cluster_results/100923_PKM2_reverse_subset_selection_train_test_loss_for_all_figs_5_repeats.csv", DataFrame)
res_df = CSV.read("Cluster_results/112723_PKM2_forward_subset_selection_train_test_loss_for_all_figs_1_repeats.csv", DataFrame)
res_df = CSV.read("Cluster_results/112723_PKM2_reverse_subset_selection_train_test_loss_for_all_figs_1_repeats.csv", DataFrame)

res_df.num_params = [15 - n_removed_param for n_removed_param in res_df.n_removed_params]

#filter res_df to remove rows with NaN in test_loss or train_loss columns
res_df_no_NaN = filter(
    row -> !isnan(row.test_loss) && !isnan(row.train_loss) && !isinf(row.train_loss) && !isinf(row.test_loss),
    res_df,
)

#drop rows of missing from res_df_no_NaN
# res_df_no_NaN = dropmissing(res_df_no_NaN)
# sort!(res_df_no_NaN, :num_params, rev = true)
# vscodedisplay(res_df_no_NaN)

#filter res_df_no_NaN to only include rows with 15 >= num_params >= 8
# res_df_no_NaN = filter(row -> row.num_params >= 8 && row.num_params <= 15, res_df_no_NaN)


##
# res_df = CSV.read("Cluster_results/112723_PKM2_forward_subset_selection_6_removed_params.csv", DataFrame)
res_df = CSV.read("Cluster_results/112723_PKM2_reverse_subset_selection_6_removed_params.csv", DataFrame)
vscodedisplay(res_df)

##
fig = Figure()
scatter_ax = Axis(
    fig[1, 1],
    # limits = (nothing, (0.00, 0.2)),
    xlabel = "Number of parameters",
    ylabel = "Test Loss",
    xticks = unique(res_df_no_NaN.num_params),
)
boxplot_ax = Axis(
    fig[1, 2],
    limits = (nothing, (-0.01, 0.6)),
    xlabel = "Number of parameters",
    ylabel = "Test Loss",
    xticks = unique(res_df_no_NaN.num_params),
)

#make dataframe with only minimal train loss for each fig, param number combination
gdf = groupby(res_df_no_NaN, :dropped_fig)
aggr_df = DataFrame()
for df in gdf
    aggr_df = vcat(
        aggr_df,
        combine(sdf -> filter(:train_loss => ≈(minimum(sdf.train_loss)), sdf), groupby(df, :num_params)),
    )
end

#plot scatter
for dropped_fig in unique(aggr_df.dropped_fig)
    println(dropped_fig)
    println(aggr_df[aggr_df.dropped_fig.==dropped_fig.&&aggr_df.num_params.==9.0, :param_subset])
    scatterlines!(
        scatter_ax,
        aggr_df[aggr_df.dropped_fig.==dropped_fig, :num_params],
        aggr_df[aggr_df.dropped_fig.==dropped_fig, :test_loss],
    )
end

#plot boxplot
mean_test = Float64[]
median_test = Float64[]
for n_params in unique(aggr_df.num_params)
    boxplot!(
        boxplot_ax,
        aggr_df[aggr_df.num_params.==n_params, :num_params],
        aggr_df[aggr_df.num_params.==n_params, :test_loss],
    )
    println(mean(aggr_df[aggr_df.num_params.==n_params, :test_loss]))
    push!(mean_test, mean(aggr_df[aggr_df.num_params.==n_params, :test_loss]))
    push!(median_test, median(aggr_df[aggr_df.num_params.==n_params, :test_loss]))
end
lines!(
    scatter_ax,
    unique(aggr_df.num_params),
    mean_test,
    color = :red,
    # linestyle = :dot,
    linewidth = 8,
    label = "Mean",
)
lines!(
    scatter_ax,
    unique(aggr_df.num_params),
    median_test,
    color = :orange,
    # linestyle = :dash,
    linewidth = 8,
    label = "Median",
)
axislegend(scatter_ax)
fig

# save("Results/$(Dates.format(now(),"mmddyy"))_subset_selection.png", fig, px_per_unit = 4)
