
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
