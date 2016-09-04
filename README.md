# ProtParam

## Installation

```
Pkg.clone("https://github.com/zmactep/ProtParam.jl.git")
```

## Usage

```
using ProtParam

rituximab = aa"QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAYPGNGDTSYNQKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWAGTTVTVSA"
protparam(rituximab)
```

## Sample output

```
ProtParam
	Query protein sequence:
		QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAYPGNGDTSYN
		QKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWAGTTVTVSA

	Theoretical pI:				8.86279296875
	Molecular weight:			12984.442140000005
	Number of amino acids in protein:	119
	Number of atoms in protein:		1781

	Total number of negatively charged residues (Asp + Glu):	7
	Total number of positively charged residues (Arg + Lys):	10

	Extinction coefficient, assuming all pairs of Cys residues form cystines:
		E:	37025
		Abs:	2.8514894672248112
	Extinction coefficient, assuming all Cys residues are reduced:
		E:	36900
		Abs:	2.841862561528576

	Half-life:
		Estimated half-life in mammalian reticulocytes, in vitro:	0.8 hour
		Estimated half-life in yeast, in vivo:				10 min
		Estimated half-life in Escherichia coli, in vivo:		>10 hour

	The instability index (II):	32.99495798319327
	This classifies the protein as stable.

	Aliphatic index of a protein:		51.68067226890756
	Grand Average of Hydropathy (GRAVY):	-0.4537815126050419

	Amino acids composition:
		A	11
		C	2
		D	4
		E	3
		F	3
		G	12
		H	1
		I	1
		K	8
		L	6
		M	3
		N	4
		P	4
		Q	7
		R	2
		T	12
		S	14
		V	8
		W	4
		Y	10

	Atoms composition:
		O	181
		H	866
		C	579
		N	150
		S	5
```
