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
    aliphatic_index

test_data = [
    (
        "rituximab",
        aa"QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAYPGNGDTSYNQKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWAGTTVTVSA",
        Dict(
            "negativecount" => 7,
            "positivecount" => 10,
            #"half_life" => ,
            "instability_index" => 32.99,
            "stability" => "stable",
            "isoelectric_point" => 8.86,
            "gravy" => -0.454,
            #"extinction_coeff" => 
            "molecular_weight" => 12984.44,
            #"absorbance" =>
            #"atom_composition" =>
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
            #"half_life" => 
            "instability_index" => 45.15,
            "stability" => "unstable",
            "isoelectric_point" => 4.57,
            "gravy" => -0.095,
            #"extinction_coeff" => 
            "molecular_weight" => 28768.78,
            #"absorbance" =>
            #"atom_composition" =>
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
            #"half_life" => 
            "instability_index" => 44.64,
            "stability" => "unstable",
            "isoelectric_point" => 10.05,
            "gravy" => -0.090,
            #"extinction_coeff" => 
            "molecular_weight" => 2250.72,
            #"absorbance" =>
            #"atom_composition" =>
            "number_of_atoms" => 322,
            "aliphatic_index" => 97.50
        )
    )
]

for (name, sequence, data) in test_data
    facts(string("Run tests for ", name)) do
        for (function_name, expected_result) in data
            result = eval(Expr(:call, symbol(function_name), sequence))
            if isa(result, Number)
                @fact result --> roughly(expected_result, 0.01)
            else
                @fact result --> expected_result
            end
        end
    end
end
