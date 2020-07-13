"Type representing atoms"
primitive type Atom <: Number 8 end

# UInt8 to Atom and back conversions
Base.convert(::Type{Atom}, atom::UInt8) = reinterpret(Atom, atom);
Base.convert(::Type{UInt8}, atom::Atom) = reinterpret(UInt8, atom);

"""
### Array of atoms
**Format**:
1. Symbol
1. Documentation string
1. Code
"""
atoms = [
    ("C",  "Carbon",   0x00),
    ("O",  "Oxigen",   0x01),
    ("H",  "Hydrogen", 0x02),
    ("N",  "Nitrogen", 0x03),
    ("S",  "Sulfur",   0x04),
];

"Prefix for all `Atom` values"
atoms_prefix = "ATOM_";

"Array of `Atom`s representation symbol, `Atom` code is a position in array"
atom_to_str = Vector{String}(undef, 0xff)

#Create all `Atom` symbols and fill their representation array
for (sym, doc, code) in atoms
    ss = Symbol(atoms_prefix, sym)
    @eval begin
        @doc $doc const $ss = convert(Atom, $code)
        atom_to_str[$code + 1] = $sym
    end
end

# Atom looks like it's symbol in String
Base.convert(::Type{String}, atom::Atom) = atom_to_str[convert(UInt8, atom) + 1]

# Atom prints as 'ATOM_<symbol>' in console
Base.show(io::IO, atom::Atom) = write(io, atoms_prefix * convert(String, atom));
