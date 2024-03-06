include("Utils.jl")

function find_all_regulator_binding_sites(
    all_regulators,
    binding_sites # list of all bindinf sites of len > 1 
)
    regulators_not_in_sites = setdiff(all_regulators, vcat(binding_sites...))
    
    all_binding_sites = [site for site in binding_sites]
    for reg in regulators_not_in_sites
    push!(all_binding_sites, [reg])
    end
    
    return all_binding_sites
end

function calculate_Z_reg(
    all_reg_binding_sites,
    kinetic_params,
    rate_data_i,
    state # active/inactive
)

    # Calculate Z^reg
    Z_reg = 1.0  # Start with 1.0 since we're multiplying
    for site in all_reg_binding_sites
        sum_term = 1.0  # The inner sum starts at 1 because of the "+1" in the formula
        for reg in site
            reg_rate = getproperty(rate_data_i, Symbol(reg)) # Get the rate data for the regulator
            K_reg = getproperty(kinetic_params, Symbol("K_", state, "_", reg, "_r"))  # Get the dissociation constant
            sum_term += reg_rate / K_reg  # Add the (regulator rate / K_regulator)
        end
        Z_reg *= sum_term  # Multiply the result by the accumulating product
    end

    return Z_reg

end

function calculate_Z_cat(S1,S2,P1, P2, K_S1,K_S2, K_P1, K_P2,
    alpha_S1_P1, 
    alpha_S1_P2, 
    alpha_S2_P1, 
    alpha_S2_P2
)
    res = (
        1 +
        (S1/K_S1) + (S2/K_S2 )+ (P1/K_P1) + (P2/K_P2)+
        ((S1*S2)/(K_S1*K_S2)) +
        ((P1*P2)/(K_P1*K_P2)) +
        (alpha_S1_P1 * (S1*P1)/(K_S1*K_P1)) +
        (alpha_S1_P2 * (S1*P2)/(K_S1*K_P2)) +
        (alpha_S2_P1 * (S2*P1)/(K_S2*K_P1)) +
        (alpha_S2_P2 * (S2*P2)/(K_S2*K_P2))
    )
    return res
end


function rate_function(
    rate_data_i,
    sub_lst,
    prod_lst,
    reg_lst,
    reg_binding_sites, 
    n, 
    kinetic_params
)

    # constant value:
    Keq = 20_000

    L = getproperty(kinetic_params, Symbol("L"))
    Vmax_a = getproperty(kinetic_params, Symbol("Vmax_a"))
    Vmax_i = getproperty(kinetic_params, Symbol("Vmax_i"))

    # regulator terms
    reg_binding_sites = find_all_regulator_binding_sites(reg_lst, reg_binding_sites)
    Z_a_reg = calculate_Z_reg(reg_binding_sites,kinetic_params,rate_data_i,"a")
    Z_i_reg = calculate_Z_reg(reg_binding_sites,kinetic_params,rate_data_i,"i")
 
    # assuming there can be only S1,S2, P1,P2
    S1 = getproperty(rate_data_i, Symbol(sub_lst[1]))
    S2 = getproperty(rate_data_i, Symbol(sub_lst[2]))
    P1 = getproperty(rate_data_i, Symbol(prod_lst[1]))
    P2 = getproperty(rate_data_i, Symbol(prod_lst[2]))

    K_a_S1 = getproperty(kinetic_params, Symbol("K_", "a", "_", sub_lst[1], "_s"))
    K_i_S1 = getproperty(kinetic_params, Symbol("K_", "i", "_", sub_lst[1], "_s"))
    K_a_S2 = getproperty(kinetic_params, Symbol("K_", "a", "_", sub_lst[2], "_s"))
    K_i_S2 = getproperty(kinetic_params, Symbol("K_", "i", "_", sub_lst[2], "_s"))
 
    K_a_P1 = getproperty(kinetic_params, Symbol("K_", "a", "_", prod_lst[1], "_p"))
    K_i_P1 = getproperty(kinetic_params, Symbol("K_", "i", "_", prod_lst[1], "_p"))
    K_a_P2 = getproperty(kinetic_params, Symbol("K_", "a", "_", prod_lst[2], "_p"))
    K_i_P2 = getproperty(kinetic_params, Symbol("K_", "i", "_", prod_lst[2], "_p"))

    alpha_S1_P1 = getproperty(kinetic_params, Symbol("alpha_", sub_lst[1], "_", prod_lst[1]))
    alpha_S1_P2 = getproperty(kinetic_params, Symbol("alpha_", sub_lst[1], "_", prod_lst[2]))
    alpha_S2_P1 = getproperty(kinetic_params, Symbol("alpha_", sub_lst[2], "_", prod_lst[1]))
    alpha_S2_P2 = getproperty(kinetic_params, Symbol("alpha_", sub_lst[2], "_", prod_lst[2]))
  
    # Calculating Z_cat terms:
    Z_a_cat = calculate_Z_cat(S1,S2,P1,P2,
        K_a_S1,K_a_S2,K_a_P1,K_a_P2,
        alpha_S1_P1, alpha_S1_P2, alpha_S2_P1, alpha_S2_P2)
    Z_i_cat = calculate_Z_cat(S1,S2,P1,P2,
        K_i_S1,K_i_S2,K_i_P1,K_i_P2,
        alpha_S1_P1, alpha_S1_P2, alpha_S2_P1, alpha_S2_P2)

    # Calculating V
    Vmax_a_rev =
        (K_a_P1 != Inf && K_a_P2 != Inf) ?
        Vmax_a * K_a_P1 * K_a_P2 / (Keq * K_a_S1 * K_a_S2) : 0.0
    Vmax_i_rev =
        (K_i_P1 != Inf && K_i_P2 != Inf) ?
        Vmax_i * K_i_P1 * K_i_P2 / (Keq * K_i_S1 * K_i_S2) : 0.0

    Rate = 
        (
            (
                Vmax_a * (S1 / K_a_S1) * (S2 / K_a_S2) -
                Vmax_a_rev * (P1 / K_a_P1) * (P2 / K_a_P2)
            ) *
            (Z_a_cat^(n-1)) *
            (Z_a_reg^n) +
            L *
            (
                Vmax_i * (S1 / K_i_S1) * (S2 / K_i_S2) -
                Vmax_i_rev * (P1 / K_i_P1) * (P2 / K_i_P2)
            ) *
            (Z_i_cat^(n-1)) *
            (Z_i_reg^n)
        ) / ((Z_a_cat^n) * (Z_a_reg^n) + L * (Z_i_cat^n) * (Z_i_reg^n))

    return Rate
end



"function to calculate rate of PKM2"
function rate_PKM2(PEP, ADP, Pyruvate, ATP, F16BP, Phenylalanine, p)
    Keq = 20_000

    Vmax_a_rev =
        (p.K_a_ATP != Inf && p.K_a_Pyruvate != Inf) ?
        p.Vmax_a * p.K_a_ATP * p.K_a_Pyruvate / (Keq * p.K_a_PEP * p.K_a_ADP) : 0.0
    Vmax_i_rev =
        (p.K_i_ATP != Inf && p.K_i_Pyruvate != Inf) ?
        p.Vmax_i * p.K_i_ATP * p.K_i_Pyruvate / (Keq * p.K_i_PEP *p. K_i_ADP) : 0.0

    Z_a_cat = (
        1 +
        (PEP / p.K_a_PEP) +
        (ATP / p.K_a_ATP) +
        (Pyruvate / p.K_a_Pyruvate) +
        (ADP / p.K_a_ADP) +
        (PEP / p.K_a_PEP) * (ADP / p.K_a_ADP) +
        (Pyruvate / p.K_a_Pyruvate) * (ATP / p.K_a_ATP) +
        p.alpha_PEP_ATP * (PEP / p.K_a_PEP) * (ATP / p.K_a_ATP) +
        p.alpha_Pyruvate_ADP * (Pyruvate / p.K_a_Pyruvate) * (ADP / p.K_a_ADP)
    )
    Z_i_cat = (
        1 +
        (PEP / p.K_i_PEP) +
        (ATP / p.K_i_ATP) +
        (Pyruvate / p.K_i_Pyruvate) +
        (ADP / p.K_i_ADP) +
        (PEP / p.K_i_PEP) * (ADP / p.K_i_ADP) +
        (Pyruvate / p.K_i_Pyruvate) * (ATP / p.K_i_ATP) +
        p.alpha_PEP_ATP * (PEP / p.K_i_PEP) * (ATP / p.K_i_ATP) +
        p.alpha_Pyruvate_ADP * (Pyruvate / p.K_i_Pyruvate) * (ADP / p.K_i_ADP)
    )
    Z_a_reg = ((1 + F16BP / p.K_a_F16BP) * (1 + Phenylalanine / p.K_a_Phenylalanine))
    Z_i_reg = ((1 + F16BP / p.K_i_F16BP) * (1 + Phenylalanine / p.K_i_Phenylalanine))

    Rate =
        (
            (
                p.Vmax_a * (PEP / p.K_a_PEP) * (ADP / p.K_a_ADP) -
                Vmax_a_rev * (Pyruvate / p.K_a_Pyruvate) * (ATP / p.K_a_ATP)
            ) *
            (Z_a_cat^3) *
            (Z_a_reg^4) +
            p.L *
            (
                p.Vmax_i * (PEP / p.K_i_PEP) * (ADP / p.K_i_ADP) -
                Vmax_i_rev * (Pyruvate / p.K_i_Pyruvate) * (ATP / p.K_i_ATP)
            ) *
            (Z_i_cat^3) *
            (Z_i_reg^4)
        ) / ((Z_a_cat^4) * (Z_a_reg^4) + p.L * (Z_i_cat^4) * (Z_i_reg^4))
    return Rate
end
