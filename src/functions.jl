"Total number of negatively charged residues (Asp + Glu)"
function negativecount(protein::AminoAcidSequence)
    comp = composition(protein)
    comp[AA_D] + comp[AA_E]
end

"Total number of positively charged residues (Arg + Lys)"
function positivecount(protein::AminoAcidSequence)
    comp = composition(protein)
    comp[AA_R] + comp[AA_K]
end

"Estimated half-life in mammalian, yeast and E. Coli"
function half_life(protein::AminoAcidSequence, make_string::Bool=true)
    result = AA_PARAMS[first(protein)].half_life
    names = Dict{ProtParamType, String}(
      PP_HL_MAM => "Mammalian, in vitro",
      PP_HL_YEAST => "Yeast, in vivo",
      PP_HL_ECOLI => "E.Coli, in vivo"
    )
    if !make_string
        result
    else
        Dict{String, String}(map(x -> names[x[1]] => x[2], collect(result)))
    end
end

"The instability index (II) of protein"
function instability_index(protein::AminoAcidSequence)
    instability = sum(DIWV((@view protein[i:i + 1])...) for i in 1:length(protein) - 1)
    round(10.0 / length(protein) * instability, digits = 2)
end

"Stability class by II"
function stability(protein::AminoAcidSequence)
    instability_index(protein) > 40 ? "unstable" : "stable"
end

"Theoretical isoelectric point of protein"
function isoelectric_point(protein::AminoAcidSequence, epsilon = 0.01)
    #these are constant pKa values
    pKa = Dict(
      AA_K => 10,
      AA_R => 12,
      AA_H => 5.98,
      AA_D => 4.05,
      AA_E => 4.45,
      AA_C => 9,
      AA_Y => 10.0
    )
    Nterm = get(Dict(
          AA_A => 7.59,
          AA_M => 7.0,
          AA_S => 6.93,
          AA_P => 8.36,
          AA_T => 6.82,
          AA_V => 7.44,
          AA_Q => 7.70
          ), first(protein), 7.5)
    Cterm = 3.55
    #"the calling function"
    charged_counts = Dict(a => count(x -> x==a, protein) for a in aa"KRHDECY")
    last_aa = last(protein)
    n = [(Cterm, 1)]
    if last_aa == AA_D
        charged_counts[last_aa] -= 1
        n = push!(n, (4.55, 1))
    elseif last_aa == AA_D
        charged_counts[last_aa] -= 1
        n = push!(n, (4.75, 1))
    end
    p = push!(map(a -> (pKa[a], charged_counts[a]), [AA_K, AA_R, AA_H]), (Nterm, 1))
    n = push!(n, map(a -> (pKa[a], charged_counts[a]), [AA_D, AA_E, AA_C, AA_Y])...)
    charge_func = pH -> begin
        CRp = [c/(1.0 + 10^(pH - pK)) for (pK, c) in p]
        CRn = [c/(1.0 + 10^(pK - pH)) for (pK, c) in n]
        sum(CRp) - sum(CRn)
    end
    current_pH = 7.0 # start from neutral pH
    current_step = 3.5
    last_charge = charge_func(current_pH)
    while abs(last_charge) >= epsilon
        current_pH =  current_pH + sign(last_charge)*current_step
        last_charge = charge_func(current_pH)
        current_step = current_step/2.0
    end
    round(current_pH, digits = 2)
end

"Grand Average of Hydropathy (GRAVY)"
function gravy(protein::AminoAcidSequence)
    round(mean([AA_PARAMS[aa].hydropathicity for aa in protein]), digits = 3)
end

"""
Pair of extinction coefficients:
1. assuming all Cys residues are reduced
1. assuming all pairs of Cys residues form cystines
"""
function extinction_coeff(protein::AminoAcidSequence)
    ext_y = 1490
    ext_w = 5500
    ext_c = 125

    comp = composition(protein)
    ec_woc = comp[AA_Y] * ext_y + comp[AA_W] * ext_w

    (ec_woc, ec_woc + div(comp[AA_C], 2) * ext_c)
end

"Molecular weight of protein"
function molecular_weight(protein::AminoAcidSequence)
    aims = [AA_MOL_PARAMS[aa].avg_isotopic_mass for aa in protein]
    round(isempty(aims) ? 0 : sum(aims) + WATER_MOL_PARAMS.avg_isotopic_mass, digits = 2)
end

"""
Absorbance 0.1% (=1 g/l)
1. assuming all Cys residues are reduced
1. assuming all pairs of Cys residues form cystines
"""
function absorbance(protein::AminoAcidSequence)
    mw = molecular_weight(protein)
    map(x -> round(x / mw, digits = 3), extinction_coeff(protein))
end

"Atoms composition"
function atom_composition(protein::AminoAcidSequence)
    atom_composition = Dict{Atom, Int}()
    for comp in [AA_MOL_PARAMS[aa].composition for aa in protein]
        for (key, val) in comp
            atom_composition[key] = get(atom_composition, key, 0) + val
        end
    end
    atom_composition[ATOM_H] = get(atom_composition, ATOM_H, 0) - 2 * (length(protein) - 1)
    atom_composition[ATOM_O] = get(atom_composition, ATOM_O, 0) - length(protein) + 1
    atom_composition
end

"Number of atoms in protein"
function number_of_atoms(protein::AminoAcidSequence)
    sum(values(atom_composition(protein)))
end

"Aliphatic index of a protein"
function aliphatic_index(protein::AminoAcidSequence)
    a = 2.9
    b = 3.9

    len = length(protein)
    comp = composition(protein)
    round(100 * (comp[AA_A] / len + a * comp[AA_V] / len + b * (comp[AA_I] + comp[AA_L]) / len), digits = 2)
end
