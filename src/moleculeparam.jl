"""
### MoleculeParam

```
MoleculeParam(carbon, hydrogen, nitrogen, oxygen, sulfur,
              monoisotopic_mass, avg_isotopic_mass)
```
**Atoms composition:**
1. `carbon`   — number of carbon atoms in molecule
1. `hydrogen` — number of hydrogen atoms in molecule
1. `nitrogen` — number of nitrogen atoms in molecule
1. `oxygen`   — number of oxygen atoms in molecule
1. `sulfur`   — number of sulfur atoms in molecule

**Weights:**
1. `monoisotopic_mass` — monoisotopic mass of molecule
1. `avg_isotopic_mass` — average isotopic mass of molecule
"""
struct MoleculeParam
    composition::Dict{Atom, Int}
    monoisotopic_mass::Float64
    avg_isotopic_mass::Float64

    function MoleculeParam(carbon::Int, hydrogen::Int, nitrogen::Int, oxygen::Int, sulfur::Int,
                           monoisotopic_mass::Float64, avg_isotopic_mass::Float64)
        new(Dict{Atom, Int}(ATOM_C => carbon,
                            ATOM_H => hydrogen,
                            ATOM_N => nitrogen,
                            ATOM_O => oxygen,
                            ATOM_S => sulfur),
            monoisotopic_mass, avg_isotopic_mass)
    end
end

"Creates `MoleculeParam` properties for all amino acids"
function create_aa_molparams()
    aa_data = Dict{AminoAcid, MoleculeParam}()
    aa_data[AA_A] = MoleculeParam(3,   7, 1, 2, 0,  71.03711,  71.0788)
    aa_data[AA_C] = MoleculeParam(3,   7, 1, 2, 1, 103.00919, 103.1388)
    aa_data[AA_D] = MoleculeParam(4,   7, 1, 4, 0, 115.02694, 115.0886)
    aa_data[AA_E] = MoleculeParam(5,   9, 1, 4, 0, 129.04259, 129.1155)
    aa_data[AA_F] = MoleculeParam(9,  11, 1, 2, 0, 147.06841, 147.1766)
    aa_data[AA_G] = MoleculeParam(2,   5, 1, 2, 0,  57.02146,  57.0519)
    aa_data[AA_H] = MoleculeParam(6,   9, 3, 2, 0, 137.05891, 137.1411)
    aa_data[AA_I] = MoleculeParam(6,  13, 1, 2, 0, 113.08406, 113.1594)
    aa_data[AA_K] = MoleculeParam(6,  14, 2, 2, 0, 128.09496, 128.1741)
    aa_data[AA_L] = MoleculeParam(6,  13, 1, 2, 0, 113.08406, 113.1594)
    aa_data[AA_M] = MoleculeParam(5,  11, 1, 2, 1, 131.04049, 131.1926)
    aa_data[AA_N] = MoleculeParam(4,   8, 2, 3, 0, 114.04293, 114.1038)
    aa_data[AA_P] = MoleculeParam(5,   9, 1, 2, 0,  97.05276,  97.1167)
    aa_data[AA_Q] = MoleculeParam(5,  10, 2, 3, 0, 128.05858, 128.1307)
    aa_data[AA_R] = MoleculeParam(6,  14, 4, 2, 0, 156.10111, 156.1875)
    aa_data[AA_S] = MoleculeParam(3,   7, 1, 3, 0,  87.03203,  87.0782)
    aa_data[AA_T] = MoleculeParam(4,   9, 1, 3, 0, 101.04768, 101.1051)
    aa_data[AA_V] = MoleculeParam(5,  11, 1, 2, 0,  99.06841,  99.1326)
    aa_data[AA_W] = MoleculeParam(11, 12, 2, 2, 0, 186.07931, 186.2132)
    aa_data[AA_Y] = MoleculeParam(9,  11, 1, 3, 0, 163.06333, 163.1760)
    aa_data
end

"`MoleculeParam` parameters for `AminoAcid`s"
AA_MOL_PARAMS = create_aa_molparams();

"`MoleculeParam` parameters for water"
WATER_MOL_PARAMS = MoleculeParam(0, 2, 0, 1, 0, 18.01056, 18.01524);
