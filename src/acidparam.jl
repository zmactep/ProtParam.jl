"""
### AminoAcidParam

```
AminoAcidParam(hl_mam, hl_yeast, hl_ecoli,
               hydropathicity, pKaC, pKaN, pKaR, pI)
```
**Half life in different conditions:**
1. `hl_mam`   — half-life in mammalian reticulocytes, in vitro
1. `hl_yeast` — half-life in yeast, in vivo
1. `hl_ecoli` — half-life in Escherichia coli, in vivo

**Hydropathicity:**
* `hydropathicity` — relative hydrophobicity or hydrophilicity of amino acid residue

**pKa:**
1. `pKaC` — pKa of amino acid C-terminus
1. `pKaN` — pKa of amino acid N-terminus
1. `pKaR` — pKa of amino acid radical
1. `pI`   — pH at which a amino acid carries no net electrical charge
"""
struct AminoAcidParam
    half_life::Dict{ProtParamType, String}
    hydropathicity::Float64
    pKaC::Float64
    pKaN::Float64
    pKaR::Float64
    pI::Float64

    function AminoAcidParam(hl_mam::String, hl_yeast::String, hl_ecoli::String,
                            hydropathicity::Float64, pKaC::Float64, pKaN::Float64, pKaR::Float64, pI::Float64)
        new(Dict{ProtParamType, String}(PP_HL_MAM => hl_mam,
                                        PP_HL_YEAST => hl_yeast,
                                        PP_HL_ECOLI => hl_ecoli),
            hydropathicity, pKaC, pKaN, pKaR, pI)
    end
end

"Creates `AminoAcidParam` properties for all amino acids"
function create_aa_params()
    aa_param = Dict{AminoAcid, AminoAcidParam}()

    aa_param[AA_A] = AminoAcidParam("4.4 hour", ">20 hour", ">10 hour",  1.8, 2.34,  9.69,   NaN,  6.00)
    aa_param[AA_R] = AminoAcidParam("1 hour",   "2 min",    "2 min",    -4.5, 2.17,  9.04, 12.48, 10.76)
    aa_param[AA_N] = AminoAcidParam("1.4 hour", "3 min",    ">10 hour", -3.5, 2.02,  8.80,   NaN,  5.41)
    aa_param[AA_D] = AminoAcidParam("1.1 hour", "3 min",    ">10 hour", -3.5, 1.88,  9.60,  3.65,  2.77)
    aa_param[AA_C] = AminoAcidParam("1.2 hour", ">20 hour", ">10 hour",  2.5, 1.96,  8.18,   NaN,  5.07)
    aa_param[AA_Q] = AminoAcidParam("0.8 hour", "10 min",   ">10 hour", -3.5, 2.17,  9.13,   NaN,  5.65)
    aa_param[AA_E] = AminoAcidParam("1 hour",   "30 min",   ">10 hour", -3.5, 2.19,  9.67,  4.25,  3.22)
    aa_param[AA_G] = AminoAcidParam("30 hour",  ">20 hour", ">10 hour", -0.4, 2.34,  9.60,   NaN,  5.97)
    aa_param[AA_H] = AminoAcidParam("3.5 hour", "10 min",   ">10 hour", -3.2, 1.82,  9.17,  6.00,  7.59)
    aa_param[AA_I] = AminoAcidParam("20 hour",  "30 min",   ">10 hour",  4.5, 2.36,  9.60,   NaN,  6.02)
    aa_param[AA_L] = AminoAcidParam("5.5 hour", "3 min",    "2 min",     3.8, 2.36,  9.60,   NaN,  5.98)
    aa_param[AA_K] = AminoAcidParam("1.3 hour", "3 min",    "2 min",    -3.9, 2.18,  8.95, 10.53,  9.74)
    aa_param[AA_M] = AminoAcidParam("30 hour",  ">20 hour", ">10 hour",  1.9, 2.28,  9.21,   NaN,  5.74)
    aa_param[AA_F] = AminoAcidParam("1.1 hour", "3 min",    "2 min",     2.8, 1.83,  9.13,   NaN,  5.48)
    aa_param[AA_P] = AminoAcidParam(">20 hour", ">20 hour", "?",        -1.6, 1.99, 10.60,   NaN,  6.30)
    aa_param[AA_S] = AminoAcidParam("1.9 hour", ">20 hour", ">10 hour", -0.8, 2.21,  9.15,   NaN,  5.68)
    aa_param[AA_T] = AminoAcidParam("7.2 hour", ">20 hour", ">10 hour", -0.7, 2.09,  9.10,   NaN,  5.60)
    aa_param[AA_W] = AminoAcidParam("2.8 hour", "3 min",    "2 min",    -0.9, 2.83,  9.39,   NaN,  5.89)
    aa_param[AA_Y] = AminoAcidParam("2.8 hour", "10 min",   "2 min",    -1.3, 2.20,  9.11,   NaN,  5.66)
    aa_param[AA_V] = AminoAcidParam("100 hour", ">20 hour", ">10 hour",  4.2, 2.32,  9.62,   NaN,  5.96)
    aa_param
end

"`AminoAcidParam` parameters for `AminoAcid`s"
AA_PARAMS = create_aa_params();
