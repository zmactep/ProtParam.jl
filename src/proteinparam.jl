"Type representing protein parameters"
bitstype 8 ProtParamType

# UInt8 to ProtParamType and back conversions
Base.convert(::Type{ProtParamType}, param::UInt8) = reinterpret(ProtParamType, param);
Base.convert(::Type{UInt8}, param::ProtParamType) = reinterpret(UInt8, param);

"""
### Array of protein parameters
**Format**:
1. Symbol
1. Documentation string
1. Code
"""
prot_params = [
    ("MW",        "Molecular weight",                                                           0x00),
    ("TOT_AA",    "Number of amino acids in protein",                                           0x01),
    ("TOT_ATOMS", "Number of atoms in protein",                                                 0x02),

    ("PI",        "Theoretical pI",                                                             0x03),

    ("COMP",      "Amino acids composition",                                                    0x04),
    ("ATOM_COMP", "Atoms composition",                                                          0x05),

    ("NEGATIVE",  "Total number of negatively charged residues (Asp + Glu)",                    0x06),
    ("POSITIVE",  "Total number of positively charged residues (Arg + Lys)",                    0x07),

    ("EXT",       "Extinction coefficient, assuming all pairs of Cys residues form cystines",    0x08),
    ("ABS",       "Absorbance 0.1% (=1 g/l), assuming all pairs of Cys residues form cystines", 0x09),
    ("EXT_NO_C",  "Extinction coefficient, assuming all Cys residues are reduced",               0x0a),
    ("ABS_NO_C",  "Absorbance 0.1% (=1 g/l), assuming all Cys residues are reduced",            0x0b),

    ("HL_MAM",    "Estimated half-life in mammalian reticulocytes, in vitro",                   0x0c),
    ("HL_YEAST",  "Estimated half-life in yeast, in vivo",                                      0x0d),
    ("HL_ECOLI",  "Estimated half-life in Escherichia coli, in vivo",                           0x0e),

    ("II",        "The instability index (II)",                                                 0x0f),
    ("II_CLASS",  "Stability class by II",                                                      0x10),

    ("ALIPH",     "Aliphatic index of a protein",                                               0x11),
    ("GRAVY",     "Grand Average of Hydropathy (GRAVY)",                                        0x12),
    ("SEQUENCE",  "Query protein sequence",                                                     0x13)
];

"Prefix for all `ProtParamType` values"
prot_params_prefix = "PP_";

"Array of `ProtParamType`s representation symbol, `ProtParamType` code is a position in array"
pp_to_str = Vector{String}(0xff)

#Create all `ProtParamType` symbols and fill their representation array
for (sym, doc, code) in prot_params
    ss = Symbol(prot_params_prefix, sym)
    @eval begin
        @doc $doc const $ss = convert(ProtParamType, $code)
        pp_to_str[$code + 1] = $sym
    end
end

# ProtParamType looks like 'PP_<symbol>' in ASCIIString
Base.convert(::Type{String}, pp::ProtParamType) = prot_params_prefix * pp_to_str[convert(UInt8, pp) + 1]

# ProtParamType prints as 'PP_<symbol>' in console
Base.show(io::IO, pp::ProtParamType) = write(io, convert(ASCIIString, pp));
