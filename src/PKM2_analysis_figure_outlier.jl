using CairoMakie, Printf
using DataFrames, CSV, StatsBase, Statistics, Dates

##
# Analyze figure Jackknife results
results_fig_loocv = CSV.read(
    "Cluster_results/102423_PKM2_fig_loocv.csv",
    DataFrame,
    #types = Float64,
)
results_fig_loocv_train_all = CSV.read(
    "Cluster_results/102523_PKM2_fig_loocv_full_train_w_all_data.csv",
    DataFrame,
    #types = Float64,
)

PKM2_data_for_fit = CSV.read("PKM2_data/PKM2_data.csv", DataFrame)
PKM2_data_for_fit.dropped_fig .= PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig

number_points_per_fig =
    combine(groupby(PKM2_data_for_fit, :dropped_fig), nrow => :number_datapoints_per_fig, renamecols = false)

min_fig_loocv = combine(sdf -> sdf[argmin(sdf.train_fit), :], groupby(results_fig_loocv, :dropped_fig))
min_fig_loocv = rename(min_fig_loocv, :test_fit => :test_loss_dropped_fig)

min_fig_loocv_train_all = combine(sdf -> sdf[argmin(sdf.train_fit), :], groupby(results_fig_loocv_train_all, :dropped_fig))
min_fig_loocv_train_all = rename(min_fig_loocv_train_all, :test_fit => :test_loss_all_figs)



analysis_df = select(innerjoin(number_points_per_fig, min_fig_loocv, min_fig_loocv_train_all; on=:dropped_fig),
    :dropped_fig,
    :number_datapoints_per_fig,
    :test_loss_all_figs,
    :test_loss_dropped_fig,
)
analysis_df.delta_test_losses = analysis_df.test_loss_dropped_fig .- analysis_df.test_loss_all_figs
analysis_df.ratio_test_losses = analysis_df.test_loss_dropped_fig ./ analysis_df.test_loss_all_figs

# vscodedisplay(analysis_df)


fig = Figure()
ax = Axis(
    fig[1, 1],
    # limits = (nothing, (0.00, 0.2)),
    xlabel = "Δ test losses w/wo figure",
    ylabel = "test loss wo figure",
)
scatter!(ax, analysis_df.delta_test_losses, analysis_df.test_loss_dropped_fig)
fig

# save("Results/$(Dates.format(now(),"mmddyy"))_outlier_figures.png", fig, px_per_unit = 4)
