module ProtParam

import Bio.Seq: AminoAcidSequence, composition,
                AminoAcid,
                AA_A, AA_C, AA_D, AA_E, AA_F,
                AA_G, AA_H, AA_I, AA_K, AA_L,
                AA_M, AA_N, AA_P, AA_Q, AA_R,
                AA_S, AA_T, AA_V, AA_W, AA_Y,
                @aa_str
import Bio.Windows: eachwindow

export Atom,
       ATOM_C, ATOM_H, ATOM_O, ATOM_N, ATOM_S,
       ProtParamHolder, protparam,
       PP_MW, PP_TOT_AA, PP_TOT_ATOMS,
       PP_PI, PP_COMP, PP_ATOM_COMP,
       PP_NEGATIVE, PP_POSITIVE,
       PP_EXT, PP_ABS, PP_EXT_NO_C, PP_ABS_NO_C,
       PP_HL_MAM, PP_HL_YEAST, PP_HL_ECOLI,
       PP_II, PP_II_CLASS,
       PP_ALIPH, PP_GRAVY,
       PP_SEQUENCE,
       @aa_str

include("atom.jl")
include("proteinparam.jl")
include("moleculeparam.jl")
include("acidparam.jl")
include("dipeptideparam.jl")
include("functions.jl")

"""
### ProtParamHolder

Holds a mapping from `ProtParamType`s to the values for the protein
"""
type ProtParamHolder
    data::Dict{ProtParamType, Any}
end

macro showdoc(x)
    :(first(first(@doc($x).content).content))
end

function Base.show(io::IO, pph::ProtParamHolder)
    write(io, "ProtParam\n")
    write(io, "\t$(@showdoc(PP_SEQUENCE)):\n")

    step = 60
    len = length(pph.data[PP_SEQUENCE])
    for i in 1:step:len
        write(io, "\t\t$(pph.data[PP_SEQUENCE][i:min(i+step - 1, len)])\n")
    end

    write(io, "\n")
    write(io, "\t$(@showdoc(PP_PI)):\t\t\t\t$(pph.data[PP_PI])\n")
    write(io, "\t$(@showdoc(PP_MW)):\t\t\t$(pph.data[PP_MW])\n")
    write(io, "\t$(@showdoc(PP_TOT_AA)):\t$(pph.data[PP_TOT_AA])\n")
    write(io, "\t$(@showdoc(PP_TOT_ATOMS)):\t\t$(pph.data[PP_TOT_ATOMS])\n")
    write(io, "\n")
    write(io, "\t$(@showdoc(PP_NEGATIVE)):\t$(pph.data[PP_NEGATIVE])\n")
    write(io, "\t$(@showdoc(PP_POSITIVE)):\t$(pph.data[PP_POSITIVE])\n")
    write(io, "\n")
    write(io, "\t$(@showdoc(PP_EXT)):\n")
    write(io, "\t\tE:\t$(pph.data[PP_EXT])\n")
    write(io, "\t\tAbs:\t$(pph.data[PP_ABS])\n")
    write(io, "\t$(@showdoc(PP_EXT_NO_C)):\n")
    write(io, "\t\tE:\t$(pph.data[PP_EXT_NO_C])\n")
    write(io, "\t\tAbs:\t$(pph.data[PP_ABS_NO_C])\n")
    write(io, "\n")
    write(io, "\tHalf-life:\n")
    write(io, "\t\t$(@showdoc(PP_HL_MAM)):\t$(pph.data[PP_HL_MAM])\n")
    write(io, "\t\t$(@showdoc(PP_HL_YEAST)):\t\t\t\t$(pph.data[PP_HL_YEAST])\n")
    write(io, "\t\t$(@showdoc(PP_HL_ECOLI)):\t\t$(pph.data[PP_HL_ECOLI])\n")
    write(io, "\n")
    write(io, "\t$(@showdoc(PP_II)):\t$(pph.data[PP_II])\n")
    write(io, "\tThis classifies the protein as $(pph.data[PP_II_CLASS]).\n")
    write(io, "\n")
    write(io, "\t$(@showdoc(PP_ALIPH)):\t\t$(pph.data[PP_ALIPH])\n")
    write(io, "\t$(@showdoc(PP_GRAVY)):\t$(pph.data[PP_GRAVY])\n")
    write(io, "\n")
    write(io, "\t$(@showdoc(PP_COMP)):\n")
    for aa in [convert(AminoAcid, c) for c in "ACDEFGHIKLMNPQRTSVWY"]
        write(io, "\t\t$(Char(aa))\t$(pph.data[PP_COMP][aa])\n")
    end
    write(io, "\n")
    write(io, "\t$(@showdoc(PP_ATOM_COMP)):\n")
    for (atom, count) in pph.data[PP_ATOM_COMP]
        write(io, "\t\t$(convert(ASCIIString, atom))\t$(count)\n")
    end
end

"Construction of protein parametes type"
function protparam(protein::AminoAcidSequence)
    data = Dict{ProtParamType, Any}()

    data[PP_COMP] = composition(protein)
    data[PP_ATOM_COMP] = atom_composition(protein)

    data[PP_MW] = molecular_weight(protein)
    data[PP_TOT_AA] = length(protein)
    data[PP_TOT_ATOMS] = number_of_atoms(protein)

    data[PP_PI] = 0

    data[PP_NEGATIVE] = negativecount(protein)
    data[PP_POSITIVE] = positivecount(protein)

    data[PP_EXT_NO_C], data[PP_EXT] = extintion_coeff(protein)
    data[PP_ABS_NO_C], data[PP_ABS] = absorbance(protein)

    hl = half_life(protein, false)
    data[PP_HL_MAM]   = hl[PP_HL_MAM]
    data[PP_HL_YEAST] = hl[PP_HL_YEAST]
    data[PP_HL_ECOLI] = hl[PP_HL_ECOLI]

    data[PP_II] = instability_index(protein)
    data[PP_II_CLASS] = stability(protein)

    data[PP_ALIPH] = aliphatic_index(protein)
    data[PP_GRAVY] = gravy(protein)

    data[PP_SEQUENCE] = protein

    ProtParamHolder(data)
end

end # module
