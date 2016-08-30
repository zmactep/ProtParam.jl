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
    names = Dict{ProtParamType, ASCIIString}(
      PP_HL_MAM => "Mammalian, in vitro",
      PP_HL_YEAST => "Yeast, in vivo",
      PP_HL_ECOLI => "E.Coli, in vivo"
    )
    if !make_string
        result
    else
        Dict{ASCIIString, Int}(map(x -> names[x[0]] => x[1], result))
    end
end

"The instability index (II) of protein"
function instability_index(protein::AminoAcidSequence)
    10.0/length(protein)*sum([DIWV(aminoacids...) for (_, aminoacids) in eachwindow(protein, 2)])
end

"Stability class by II"
function stability(protein::AminoAcidSequence)
    instability_index(protein) > 40 ? "unstable" : "stable"
end

"Theoretical isoelectric point of protein"
function isoelectric_point(protein::AminoAcidSequence)
    0
end

"Grand Average of Hydropathy (GRAVY)"
function gravy(protein::AminoAcidSequence)
    mean([AA_PARAMS[aa].hydropathicity for aa in protein])
end

"""
Pair of extintion coefficients:
1. assuming all Cys residues are reduced
1. assuming all pairs of Cys residues form cystines
"""
function extintion_coeff(protein::AminoAcidSequence)
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
    isempty(aims) ? 0 : sum(aims) + WATER_MOL_PARAMS.avg_isotopic_mass
end

"""
Absorbance 0.1% (=1 g/l)
1. assuming all Cys residues are reduced
1. assuming all pairs of Cys residues form cystines
"""
function absorbance(protein::AminoAcidSequence)
    mw = molecular_weight(protein)
    map(x -> x / mw, extintion_coeff(protein))
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
    100 * (comp[AA_A] / len + a * comp[AA_V] / len + b * (comp[AA_I] + comp[AA_L]) / len)
end
