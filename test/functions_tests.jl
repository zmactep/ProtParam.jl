# This test set ensures that methods from `src/functions.jl` behave the same way as they should
# This is checked by comparing results of particular methods with corresponding output from ExPASy.
import ProtParam: *
rituximab = aa"QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAYPGNGDTSYNQKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWAGTTVTVSA"
PCNA = aa"MFEARLVQGSILKKVLEALKDLINEACWDISSSGVNLQSMDSSHVSLVQLTLRSEGFDTYRCDRNLAMGVNLTSMSKILKCAGNEDIITLRAEDNADTLALVFEAPNQEKVSDYEMKLMDLDVEQLGIPEQEYSCVVKMPSGEFARICRDLSHIGDAVVISCAKDGVKFSASGELGNGNIKLSQTSNVDKEEEAVTIEMNEPVQLTFALRYLNFFTKATPLSSTVTLSMSADVPLVVEYKIADMGHLKYYLAPKIEDEEGS" 

facts("negativecount computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact negativecount(protein) --> roughly(expected_result, 0.01)
    end
end

facts("positivecount computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact positivecount(protein) --> roughly(expected_result, 0.01)
    end
end

facts("half_life computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact half_life(protein) --> roughly(expected_result, 0.01)
        @pending half_life(protein, false)
    end
end

facts("instability_index computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact instability_index(protein) --> roughly(expected_result, 0.01)
    end
end

facts("stability computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact stability(protein) --> roughly(expected_result, 0.01)
    end
end

# compute theoretical pI for several sequences and compare with Compute pI/Mw output
facts("isoelectric_point computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact isoelectric_point(protein) --> roughly(expected_result, 0.01)
    end
end

facts("gravy computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact gravy(protein) --> roughly(expected_result, 0.01)
    end
end

facts("extinction_coeff computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact extinction_coeff(protein) --> roughly(expected_result, 0.01)
    end
end

facts("molecular_weight computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact molecular_weight(protein) --> roughly(expected_result, 0.01)
    end
end

facts("absorbance computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact absorbance(protein) --> roughly(expected_result, 0.01)
    end
end

facts("atom_composition computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact atom_composition(protein) --> roughly(expected_result, 0.01)
    end
end

facts("number_of_atoms computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact number_of_atoms(protein) --> roughly(expected_result, 0.01)
    end
end

facts("aliphatic_index computation") do
    for (protein, expected_result) in [
        (aa"RCDRNLAMGVNLTSMSKILK", 10.05),
        (rituximab, 8.86),
        (PCNA, 4.57)
    ]
        @fact aliphatic_index(protein) --> roughly(expected_result, 0.01)
    end
end
