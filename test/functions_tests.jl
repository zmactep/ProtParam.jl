# This test set ensures that methods from `src/functions.jl` behave the same way as they should
# This is checked by comparing results of particular methods with corresponding output from ExPASy.
import ProtParam: 
    negativecount,
    positivecount,
    half_life,
    instability_index,
    stability,
    isoelectric_point, 
    gravy,
    extinction_coeff,
    molecular_weight,
    absorbance,
    atom_composition,
    number_of_atoms,
    aliphatic_index,
    ProtParamType

test_data = [
    (
        "rituximab",
        aa"QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAYPGNGDTSYNQKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWAGTTVTVSA",
        Dict(
            "negativecount" => 7,
            "positivecount" => 10,
            "half_life" => Dict(
                "Mammalian, in vitro" => "0.8 hour",
                "Yeast, in vivo" => "10 min",
                "E.Coli, in vivo" => "10 hour"
            ), # this dict is used to test half_life function with `make_string` parameter set to true
            "instability_index" => 32.99,
            "stability" => "stable",
            "isoelectric_point" => 8.86,
            "gravy" => -0.454,
            "extinction_coeff" => (36900,37025), 
            "molecular_weight" => 12984.44,
            "absorbance" => (2.842, 2.851),
            "atom_composition" => Dict(
                ATOM_S => 5,
                ATOM_C => 579,
                ATOM_H => 866,
                ATOM_N => 150,
                ATOM_O => 181
            ),
            "number_of_atoms" => 1781,
            "aliphatic_index" => 51.68
        )
    ), 
    (   
        "PCNA",
        aa"MFEARLVQGSILKKVLEALKDLINEACWDISSSGVNLQSMDSSHVSLVQLTLRSEGFDTYRCDRNLAMGVNLTSMSKILKCAGNEDIITLRAEDNADTLALVFEAPNQEKVSDYEMKLMDLDVEQLGIPEQEYSCVVKMPSGEFARICRDLSHIGDAVVISCAKDGVKFSASGELGNGNIKLSQTSNVDKEEEAVTIEMNEPVQLTFALRYLNFFTKATPLSSTVTLSMSADVPLVVEYKIADMGHLKYYLAPKIEDEEGS",
        Dict(
            "negativecount" => 41,
            "positivecount" => 24,
            "half_life" => Dict(
                "Mammalian, in vitro" => "30 hour",
                "Yeast, in vivo" => ">20 hour",
                "E.Coli, in vivo" => ">10 hour"
            ), # this dict is used to test half_life function with `make_string` parameter set to true
            "instability_index" => 45.15,
            "stability" => "unstable",
            "isoelectric_point" => 4.57,
            "gravy" => -0.095,
            "extinction_coeff" => (15930,16305), 
            "molecular_weight" => 28768.78,
            "absorbance" => (0.554, 0.567),
            "atom_composition" => Dict(
                ATOM_S => 16,
                ATOM_C => 1257,
                ATOM_H => 2020,
                ATOM_N => 328,
                ATOM_O => 408
            ),
            "number_of_atoms" => 4029,
            "aliphatic_index" => 94.87
        )
    ),
    (
        "PCNA fragment",
        aa"RCDRNLAMGVNLTSMSKILK",
        Dict(
            "negativecount" => 1,
            "positivecount" => 4,
            "half_life" => Dict(
                "Mammalian, in vitro" => "1 hour",
                "Yeast, in vivo" => "2 min",
                "E.Coli, in vivo" => "2 min"
            ), # this dict is used to test half_life function with `make_string` parameter set to true
            "instability_index" => 44.64,
            "stability" => "unstable",
            "isoelectric_point" => 10.05,
            "gravy" => -0.090,
            "extinction_coeff" => (0, 0),
            "molecular_weight" => 2250.72,
            "absorbance" => (0, 0),
            "atom_composition" => Dict(
                ATOM_S => 3,
                ATOM_C => 93,
                ATOM_H => 168,
                ATOM_N => 30,
                ATOM_O => 28
            ),
            "number_of_atoms" => 322,
            "aliphatic_index" => 97.50
        )
    )
]

aa1 = test_data[1][2]
println(atom_composition(aa1))
#println(show(half_life(aa1)))
for (name, sequence, data) in test_data
    facts(string("Run tests for ", name)) do
        for (function_name, expected_result) in data
            result = eval(Expr(:call, symbol(function_name), sequence))
            if isa(result, Number)
                @fact result --> roughly(expected_result, 0.01)
            elseif isa(result, Tuple)
                for (u, v) in zip(result, expected_result)
                    @fact u --> (isa(v, Number) ? roughly(v, 0.01) : v)
                end
            else
                @fact result --> expected_result
            end
        end
    end
end
