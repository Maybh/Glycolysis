using CairoMakie, CSV, DataFrames, Dates, Statistics, Printf, StatsBase, Distributions

#using FileIO
#CairoMakie.activate!(type = "png")


##
#=
Describe bootstrap results statistics
=#

# Load and process bootstrap data
bootstrap_data = CSV.read(
    # "Cluster_results/080723_PKM2_fitting_10000_fig_bootstrap_w_5_repeats_rand_x0_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf.csv",
    # "Cluster_results/100923_PKM2_fitting_10000_fig_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_1_alphaPyrADP_1.csv",
    # "Cluster_results/100923_PKM2_fitting_10000_fig_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_0_alphaPyrADP_0.csv",
    "Cluster_results/100923_PKM2_fitting_10000_fig_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_0_alphaPyrADP_1.csv",
    DataFrame,
)
# Only include data with minimal fitness
bootstrap_data = combine(
    sdf -> filter(:fitness => ≈(minimum(sdf.fitness)), sdf),
    groupby(bootstrap_data, :bootrap_fig_list),
)


# Only include data that has at least one figure with [Pyruvate] != 0.0
PKM2_data = CSV.read("PKM2_data/PKM2_data.csv", DataFrame)
PKM2_data.source .= PKM2_data.Article .* "_" .* PKM2_data.Fig
Pyruvate_figs = unique(PKM2_data[PKM2_data.Pyruvate.>0.0, :].source)
bootstrap_data = subset(bootstrap_data, :bootrap_fig_list => ByRow(x -> any(occursin.(Pyruvate_figs, x))))

# Only include data that has at least one figure with [F16BP] < 10 µM
# F16BP_figs = unique(PKM2_data[PKM2_data.F16BP .> 0.0 .&& PKM2_data.F16BP .< 1e-5, :].source)
# bootstrap_data = subset(bootstrap_data, :bootrap_fig_list => ByRow(x -> any(occursin.(F16BP_figs, x))))

# Select only the parameters that are not fixed by Ka=Ki or K=Inf
params_for_plot = [:L, :K_a_PEP, :K_i_PEP, :K_a_ADP, :K_a_Pyruvate, :K_a_ATP, :K_a_F16BP, :K_i_Phenylalanine]
data_for_describe = bootstrap_data[:, params_for_plot]

# Calculate mean, std, and CV for each parameter assuming Log Normal or Normal distribution
describe(
    data_for_describe,
    (x -> mean(fit_mle(Normal, x))) => :mean,
    (x -> std(fit_mle(Normal, x))) => :std,
    (x -> std(fit_mle(Normal, x)) / mean(fit_mle(Normal, x))) => :CV,
)

vscodedisplay(describe(
    data_for_describe,
    (x -> mean(fit_mle(LogNormal, x))) => :mean,
    (x -> std(fit_mle(LogNormal, x))) => :std,
    (x -> std(fit_mle(LogNormal, x)) / mean(fit_mle(LogNormal, x))) => :CV,
))

##
# Plot heatmap of correlation between model parameters

# Load and process data with repeats
data_for_plot = CSV.read(
    # "Cluster_results/080723_PKM2_fitting_10000_fig_bootstrap_w_5_repeats_rand_x0_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf.csv",
    # "Cluster_results/100923_PKM2_fitting_10000_fig_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_1_alphaPyrADP_1.csv",
    "Cluster_results/100923_PKM2_fitting_10000_fig_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_0_alphaPyrADP_1.csv",
    DataFrame,
)
# Only include data with minimal fitness
data_for_plot = combine(
    sdf -> filter(:fitness => ≈(minimum(sdf.fitness)), sdf),
    groupby(data_for_plot, :bootrap_fig_list),
)
# Only include data that has at least one figure with [Pyruvate] != 0.0
Pyruvate_figs = unique(PKM2_data[PKM2_data.Pyruvate.>0.0, :].source)
data_for_plot = subset(data_for_plot, :bootrap_fig_list => ByRow(x -> any(occursin.(Pyruvate_figs, x))))

# Only include data that has at least one figure with [F16BP] < 10 µM
# F16BP_figs = unique(PKM2_data[PKM2_data.F16BP .> 0.0 .&& PKM2_data.F16BP .< 1e-5, :].source)
# data_for_plot = subset(data_for_plot, :bootrap_fig_list => ByRow(x -> any(occursin.(F16BP_figs, x))))

params_for_plot = [:L, :K_a_PEP, :K_a_ADP, :K_a_Pyruvate, :K_a_ATP, :K_a_F16BP, :K_i_PEP, :K_i_Phenylalanine]
cov(Matrix(data_for_plot[:, params_for_plot]))
function plot_param_corr_heatmap(data_for_plot, fig = Figure())
    ax, hm = heatmap(
        fig[1, 1],
        cor(log.(Matrix(data_for_plot))),
        colormap = :balance,
        colorrange = (-1, 1),
        axis = (
            xticks = (1:length(names(data_for_plot)), names(data_for_plot)),
            yticks = (1:length(names(data_for_plot)), names(data_for_plot)),
            xticklabelrotation = π / 2,
        ),
    )
    Colorbar(
        fig[1, 2],
        hm,
        height = Relative(2 / 3),
        label = "Pearson correlation",
        flip_vertical_label = false,
    )
    fig
end
# vscodedisplay(cor(Matrix(bootstrap_data[:, params_for_plot])))
plot_param_corr_heatmap(data_for_plot[:, params_for_plot])






##
# Plot histogram bootstrapping results

# Load and process data with repeats
data_for_plot = CSV.read(
    # "Cluster_results/080723_PKM2_fitting_10000_fig_bootstrap_w_5_repeats_rand_x0_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf.csv",
    # "Cluster_results/100923_PKM2_fitting_10000_fig_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_1_alphaPyrADP_1.csv",
    "Cluster_results/100923_PKM2_fitting_10000_fig_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_0_alphaPyrADP_1.csv",
    DataFrame,
)
# Only include data with minimal fitness
data_for_plot = combine(
    sdf -> filter(:fitness => ≈(minimum(sdf.fitness)), sdf),
    groupby(data_for_plot, :bootrap_fig_list),
)
# Only include data that has at least one figure with [Pyruvate] != 0.0
Pyruvate_figs = unique(PKM2_data[PKM2_data.Pyruvate.>0.0, :].source)
data_for_plot = subset(data_for_plot, :bootrap_fig_list => ByRow(x -> any(occursin.(Pyruvate_figs, x))))

# Only include data that has at least one figure with [F16BP] < 10 µM
# F16BP_figs = unique(PKM2_data[PKM2_data.F16BP .> 0.0 .&& PKM2_data.F16BP .< 1e-5, :].source)
# data_for_plot = subset(data_for_plot, :bootrap_fig_list => ByRow(x -> any(occursin.(F16BP_figs, x))))
params_for_plot = [:L, :K_a_PEP, :K_i_PEP, :K_a_ADP, :K_a_Pyruvate, :K_a_ATP, :K_a_F16BP, :K_i_Phenylalanine]
bootstrap_data = data_for_plot[:, params_for_plot]

# # Generate and sample params from Multivariate Log Normal Distribution MLE fits of data to confirm that it is a good distribution to use
# params_for_plot = [:L, :K_a_PEP, :K_i_PEP, :K_a_ADP, :K_a_Pyruvate, :K_a_ATP, :K_a_F16BP, :K_i_Phenylalanine]
# data_for_describe = data_for_plot[:, params_for_plot]
# μ = vec(mean(log.(Matrix(data_for_describe[:, params_for_plot])), dims = 1))
# Σ = cov(log.(Matrix(data_for_describe[:, params_for_plot])))
# # Create a bivariate log-normal distribution with the specified means and covariance
# dist = MvLogNormal(μ, Σ)

# # Generate 1000 samples from the dist distribution and put into DataFrame with params_for_plot columns
# bootstrap_data = DataFrame(rand(dist, 10_000)', params_for_plot)

function plot_param_bootsrap_histograms(data_for_plot; fig = Figure(), num_columns = 3)
    axs = []
    for (i, fit_parameter) in enumerate(names(data_for_plot))
        hist(
            fig[
                i % num_columns == 0 ? i ÷ num_columns : 1 + i ÷ num_columns,
                i % num_columns == 0 ? num_columns : i % num_columns,
            ],
            data_for_plot[!, fit_parameter],
            # setting for all params except weights
            bins = begin
                if any(contains.(fit_parameter, ["K_"]))
                    exp10.(range(-8, stop = 3, length = 100))
                elseif any(contains.(fit_parameter, ["L"]))
                    exp10.(range(-4, stop = 4, length = 100))
                elseif any(contains.(fit_parameter, ["V"]))
                    exp10.(range(-3, stop = 3, length = 100))
                end
            end,
            axis = (
                xlabel = replace(fit_parameter, "_" => " ") * ", M",
                xscale = log10,
                xlabelpadding = 1,
                xticklabelpad = 0,
                yticklabelpad = 0,
                #xticks = LogTicks(Makie.IntegerTicks()),
                #xticks = WilkinsonTicks(4)
            ),
            normalization = :probability,
        )
        push!(axs, current_axis())
    end
    #linkyaxes!(axs...)
    #[hideydecorations!(ax, grid = false) for ax in axs[[2:4; 6:end]]]
    [hideydecorations!(ax) for ax in axs]
    fig
end
# plot_param_bootsrap_histograms(bootstrap_data; num_columns=4, fig=Figure(resolution=(2000, 2000)))

fig = Figure()
plot_param_bootsrap_histograms(bootstrap_data; num_columns = 4, fig = fig[1, 1])
fig





##
# include various functions used for fitting
# Make sure rate_PKM2 function is consistent with the one used for fitting by setting relevant Ka=Ki or K=Inf
include("PKM2_fitting_functions.jl")

bootstrap_param_data = CSV.read("Results/100923_PKM2_fitting_8_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf.csv", DataFrame)
bootstrap_param_data = CSV.read("Results/100923_PKM2_fitting_8_bootstrap_a_neq_i_PEP_Phe_F16BP_aPhe_iF16BP_Inf_alphaPEPATP_0.csv", DataFrame)
bootstrap_param_data = CSV.read("Results/102423_PKM2_fitting_4_bootstrap_full_model.csv", DataFrame)

raw_data = CSV.read("PKM2_data/PKM2_data.csv", DataFrame)
raw_data.source .= raw_data.Article .* "_" .* raw_data.Fig

#Filter data. Make sure its consistent with data used for fitting.
# filter!((row -> !(row.source == "Dombrauckas_Biochem_2005_Fig5A" && (row.PEP < 0.1e-3 || row.Rate < 0.2))), raw_data)

function unit_conc_convertion(conc::Number)
    if conc >= 1e-2
        str = @sprintf("%.0f", conc / 1e-3) * "mM"
    elseif conc >= 1e-3
        str = @sprintf("%.1f", conc / 1e-3) * "mM"
    elseif conc >= 1e-5
        str = @sprintf("%.0f", conc / 1e-6) * "µM"
    elseif conc >= 1e-6
        str = @sprintf("%.1f", conc / 1e-6) * "µM"
    elseif conc >= 1e-8
        str = @sprintf("%.0f", conc / 1e-9) * "nM"
    elseif conc >= 1e-9
        str = @sprintf("%.1f", conc / 1e-9) * "nM"
    elseif conc == 0.0
        str = "0.0     "
    else
        str = "Number Out of Scale"
    end
    return str
end

"function for plotting fit of `rate_PKM2()` with kinetic parameters of fit from `bootstrap_param_data` on `raw_data`"
function plot_fit_on_data(
    bootstrap_param_data,
    raw_data,
    param_names;
    fig = Figure(),
    num_col = 5,
    legend_font_size = 12,
    markersize = 3,
    linewidth = 1,
)
    kinetic_params = Vector(
        eachrow(bootstrap_param_data[:, param_names][bootstrap_param_data.fitness.==minimum(bootstrap_param_data.fitness), :])[1],
    )
    # kinetic_params = [mean(Matrix(bootstrap_param_data[:, Between(:L, names(bootstrap_param_data)[end])]), dims = 1)...]
    
    # Process raw_data to calculate figure weights and correct the rate data
    # Add a new column to data to assign an integer to each source/figure from publication
    raw_data.fig_num = vcat(
        [
            i * ones(Int64, count(==(unique(raw_data.source)[i]), raw_data.source)) for
            i = 1:length(unique(raw_data.source))
        ]...,
    )

    # Make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(raw_data.fig_num .== i) for i in unique(raw_data.fig_num)]
    figure_weights = figure_weights_PKM2(kinetic_params, raw_data, fig_point_indexes)
    raw_data.Rate .= raw_data.Rate ./ figure_weights

    axs = []
    for (i, Fig) in enumerate(unique(raw_data.source))
        PEP = raw_data[raw_data.source.==Fig, :PEP]
        ADP = raw_data[raw_data.source.==Fig, :ADP]
        Pyruvate = raw_data[raw_data.source.==Fig, :Pyruvate]
        ATP = raw_data[raw_data.source.==Fig, :ATP]
        F16BP = raw_data[raw_data.source.==Fig, :F16BP]
        Phenylalanine = raw_data[raw_data.source.==Fig, :Phenylalanine]
        Rate = raw_data[raw_data.source.==Fig, :Rate]
        X_axis_label = raw_data[raw_data.source.==Fig, :X_axis_label]
        
        #make sure the X_axis_label is the same for all points in the figure
        @assert all(X_axis_label .== X_axis_label[1])
        
        gd =
            fig[i % num_col == 0 ? i ÷ num_col : 1 + i ÷ num_col, i % num_col == 0 ? num_col : i % num_col] =
                GridLayout()
        ax = Axis(
            gd[1, 1],
            xticks = LinearTicks(2),
            yticks = LinearTicks(2),
            yticklabelrotation = π / 2,
            xlabelpadding = 1,
            xticklabelpad = 0,
            yticklabelpad = 0,
        )
        if i % num_col == 1
            ax.ylabel = "PKM2 Rate"
        end
        scas = []
        lins = []
        sca =nothing
        lin=nothing
        concs = []
        if X_axis_label[1] == "F16BP"
                ax.xlabel = "F16BP, M"
            for unique_PEP in unique(PEP),
                unique_ATP in unique(ATP),
                unique_Pyruvate in unique(Pyruvate),
                unique_ADP in unique(ADP),
                phenylalanine in unique(Phenylalanine)

                if length(
                    Rate[(PEP.==unique_PEP).&(ADP.==unique_ADP).&(ATP.==unique_ATP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                ) .> 0
                    global sca = scatter!(
                        F16BP[(PEP.==unique_PEP).&(ADP.==unique_ADP).&(ATP.==unique_ATP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        Rate[(PEP.==unique_PEP).&(ADP.==unique_ADP).&(ATP.==unique_ATP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        markersize = markersize,
                    )
                    global lin = lines!(
                        0 .. maximum(F16BP),
                        x -> rate_PKM2(unique_PEP, unique_ADP, unique_Pyruvate, unique_ATP, x, phenylalanine, kinetic_params),
                        linewidth = linewidth,
                    )
                end
                concs = vcat(
                    concs,
                    (; PEP = unique_PEP, ATP = unique_ATP, Pyruvate = unique_Pyruvate, ADP = unique_ADP, Phenylalanine = phenylalanine),
                )
                push!(scas, sca)
                push!(lins, lin)
            end
        elseif X_axis_label[1] == "ADP"
                ax.xlabel = "ADP, M"
            for unique_F16BP in unique(F16BP),
                unique_PEP in unique(PEP),
                unique_ATP in unique(ATP),
                unique_Pyruvate in unique(Pyruvate),
                phenylalanine in unique(Phenylalanine)

                if length(
                    Rate[(PEP.==unique_PEP).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                ) .> 0
                    global sca = scatter!(
                        ADP[(PEP.==unique_PEP).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        Rate[(PEP.==unique_PEP).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        markersize = markersize,
                    )
                    global lin = lines!(
                        0 .. maximum(ADP),
                        x -> rate_PKM2(unique_PEP, x, unique_Pyruvate, unique_ATP, unique_F16BP, phenylalanine, kinetic_params),
                        linewidth = linewidth,
                    )
                    concs = vcat(
                        concs,
                        (;
                            PEP = unique_PEP,
                            ATP = unique_ATP,
                            Pyruvate = unique_Pyruvate,
                            F16BP = unique_F16BP,
                            Phenylalanine = phenylalanine,
                        ),
                    )
                    push!(scas, sca)
                    push!(lins, lin)
                end
            end
        elseif X_axis_label[1] == "PEP"
                ax.xlabel = "PEP, M"
            for unique_F16BP in unique(F16BP),
                unique_ADP in unique(ADP),
                unique_ATP in unique(ATP),
                unique_Pyruvate in unique(Pyruvate),
                phenylalanine in unique(Phenylalanine)

                if length(
                    Rate[(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                ) .> 0
                    global sca = scatter!(
                        PEP[(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        Rate[(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        markersize = markersize,
                    )
                    global lin = lines!(
                        0 .. maximum(PEP),
                        x -> rate_PKM2(x, unique_ADP, unique_Pyruvate, unique_ATP, unique_F16BP, phenylalanine, kinetic_params),
                        linewidth = linewidth,
                    )
                    concs = vcat(
                        concs,
                        (;
                            ATP = unique_ATP,
                            Pyruvate = unique_Pyruvate,
                            ADP = unique_ADP,
                            F16BP = unique_F16BP,
                            Phenylalanine = phenylalanine,
                        ),
                    )
                    push!(scas, sca)
                    push!(lins, lin)
                end
            end
        elseif X_axis_label[1] == "ATP"
                ax.xlabel = "ATP, M"
            for unique_F16BP in unique(F16BP),
                unique_ADP in unique(ADP),
                unique_PEP in unique(PEP),
                unique_Pyruvate in unique(Pyruvate),
                phenylalanine in unique(Phenylalanine)

                if length(
                    Rate[(F16BP.==unique_F16BP).&(PEP.==unique_PEP).&(ADP.==unique_ADP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                ) .> 0
                    global sca = scatter!(
                        ATP[(F16BP.==unique_F16BP).&(PEP.==unique_PEP).&(ADP.==unique_ADP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        Rate[(F16BP.==unique_F16BP).&(PEP.==unique_PEP).&(ADP.==unique_ADP).&(Pyruvate.==unique_Pyruvate).&(Phenylalanine.==phenylalanine)],
                        markersize = markersize,
                    )
                    global lin = lines!(
                        0 .. maximum(ATP),
                        x -> rate_PKM2(unique_PEP, unique_ADP, unique_Pyruvate, x, unique_F16BP, phenylalanine, kinetic_params),
                        linewidth = linewidth,
                    )
                    concs = vcat(
                        concs,
                        (;
                            PEP = unique_PEP,
                            Pyruvate = unique_Pyruvate,
                            ADP = unique_ADP,
                            F16BP = unique_F16BP,
                            Phenylalanine = phenylalanine,
                        ),
                    )
                    push!(scas, sca)
                    push!(lins, lin)
                end
            end
        elseif X_axis_label[1] == "Pyruvate"
                ax.xlabel = "Pyruvate, M"
            for unique_F16BP in unique(F16BP),
                unique_PEP in unique(PEP),
                unique_ATP in unique(ATP),
                unique_ADP in unique(ADP),
                phenylalanine in unique(Phenylalanine)

                if length(
                    Rate[(PEP.==unique_PEP).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP).&(Phenylalanine.==phenylalanine)],
                ) .> 0
                    global sca = scatter!(
                        Pyruvate[(PEP.==unique_PEP).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP).&(Phenylalanine.==phenylalanine)],
                        Rate[(PEP.==unique_PEP).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP).&(Phenylalanine.==phenylalanine)],
                        markersize = markersize,
                    )
                    global lin = lines!(
                        0 .. maximum(Pyruvate),
                        x -> rate_PKM2(unique_PEP, unique_ADP, x, unique_ATP, unique_F16BP, phenylalanine, kinetic_params),
                        linewidth = linewidth,
                    )
                    concs = vcat(
                        concs,
                        (; PEP = unique_PEP, ATP = unique_ATP, ADP = unique_ADP, F16BP = unique_F16BP, Phenylalanine = phenylalanine),
                    )
                    push!(scas, sca)
                    push!(lins, lin)
                end
            end
        elseif X_axis_label[1] == "Phenylalanine"
            ax.xlabel = "Phenylalanine, M"
            for unique_PEP in unique(PEP),
                unique_ATP in unique(ATP),
                unique_Pyruvate in unique(Pyruvate),
                unique_ADP in unique(ADP),
                unique_F16BP in unique(F16BP)

                if length(
                    Rate[(PEP.==unique_PEP).&(Pyruvate.==unique_Pyruvate).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP)],
                ) .> 0
                    global sca = scatter!(
                        Phenylalanine[(PEP.==unique_PEP).&(Pyruvate.==unique_Pyruvate).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP)],
                        Rate[(PEP.==unique_PEP).&(Pyruvate.==unique_Pyruvate).&(F16BP.==unique_F16BP).&(ATP.==unique_ATP).&(ADP.==unique_ADP)],
                        markersize = markersize,
                    )
                    global lin = lines!(
                        0 .. maximum(Phenylalanine),
                        x -> rate_PKM2(unique_PEP, unique_ADP, unique_Pyruvate, unique_ATP, unique_F16BP, x, kinetic_params),
                        linewidth = linewidth,
                    )
                    concs =
                        vcat(concs, (; PEP = unique_PEP, ATP = unique_ATP, Pyruvate = unique_Pyruvate, ADP = unique_ADP, F16BP = unique_F16BP))
                    push!(scas, sca)
                    push!(lins, lin)
                end
            end
        else
            Axis(
                fig[
                    i % num_col == 0 ? i ÷ num_col : 1 + i ÷ num_col,
                    i % num_col == 0 ? num_col : i % num_col,
                ],
                xlabel = "Something went wrong and this Fig didn't plot",
                #title = "$Fig",
            )
        end
        Legend(
            gd[1, 2],
            lins,
            begin
                strs = Vector{String}()
                for i = 1:length(lins)
                    str = ""
                    if length(unique(PEP)) > 1 && ax.xlabel.val != "PEP, M"
                        str = str * unit_conc_convertion(concs[i].PEP) * "    "
                    end
                    if length(unique(ATP)) > 1 && ax.xlabel.val != "ATP, M"
                        str = str * unit_conc_convertion(concs[i].ATP) * "    "
                    end
                    if length(unique(Pyruvate)) > 1 && ax.xlabel.val != "Pyruvate, M"
                        str = str * unit_conc_convertion(concs[i].Pyruvate) * "    "
                    end
                    if length(unique(ADP)) > 1 && ax.xlabel.val != "ADP, M"
                        str = str * unit_conc_convertion(concs[i].ADP) * "    "
                    end
                    if length(unique(F16BP)) > 1 && ax.xlabel.val != "F16BP, M"
                        str = str * unit_conc_convertion(concs[i].F16BP) * "    "
                    end
                    if length(unique(Phenylalanine)) > 1 && ax.xlabel.val != "Phenylalanine, M"
                        str = str * unit_conc_convertion(concs[i].Phenylalanine) * "    "
                    end
                    push!(strs, str)
                end
                strs
            end,
            begin
                str1 = ""
                str2 = "      "
                if length(unique(PEP)) == 1
                    PEP[1] != 0.0 && (str1 = str1 * "PEP = " * unit_conc_convertion(unique(PEP)[1]) * "\n")
                elseif ax.xlabel.val != "PEP, M"
                    str2 = str2 * "PEP" * "       "
                end
                if length(unique(ATP)) == 1
                    ATP[1] != 0.0 && (str1 = str1 * "ATP = " * unit_conc_convertion(unique(ATP)[1]) * "\n")
                elseif ax.xlabel.val != "ATP, M"
                    str2 = str2 * "ATP" * "       "
                end
                if length(unique(Pyruvate)) == 1
                    Pyruvate[1] != 0.0 &&
                        (str1 = str1 * "Pyruvate = " * unit_conc_convertion(unique(Pyruvate)[1]) * "\n")
                elseif ax.xlabel.val != "Pyruvate, M"
                    str2 = str2 * "Pyruvate" * "       "
                end
                if length(unique(ADP)) == 1
                    ADP[1] != 0.0 && (str1 = str1 * "ADP = " * unit_conc_convertion(unique(ADP)[1]) * "\n")
                elseif ax.xlabel.val != "ADP, M"
                    str2 = str2 * "ADP" * "       "
                end
                if length(unique(F16BP)) == 1
                    F16BP[1] != 0.0 &&
                        (str1 = str1 * "F16BP = " * unit_conc_convertion(unique(F16BP)[1]) * "\n")
                elseif ax.xlabel.val != "F16BP, M"
                    str2 = str2 * "F16BP" * "       "
                end
                if length(unique(Phenylalanine)) == 1
                    Phenylalanine[1] != 0.0 &&
                        (str1 = str1 * "Phe = " * unit_conc_convertion(unique(Phenylalanine)[1]) * "\n")
                elseif ax.xlabel.val != "Phenylalanine, M"
                    str2 = str2 * "Phe" * "       "
                end
                str1 * "\n" * str2
            end;
            patchsize = (2.0f0, 2.0f0),
            patchlabelgap = 2,
            padding = (1.0f0, 1.0f0, 1.0f0, 1.0f0),
            framevisible = false,
            rowgap = 0,
            titlegap = 2,
            titlehalign = :left,
            titlesize = legend_font_size,
            labelsize = legend_font_size,
        )
        push!(axs, current_axis())
        colgap!(gd, 1)
        Label(gd[1, :, Top()], replace(Fig, "_" => " "), padding = (0, 0, 3, 0))
    end
    fig
end


plot_fit_on_data(
    bootstrap_param_data,
    raw_data,
    param_names;
    fig = Figure(resolution = (3000 / 1.5, 2000 / 1.5)),
    num_col = 6,
)


