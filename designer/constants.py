import pandas as pd

from designer.paths import DATA_DIR


CANONICAL_AMINO_ACIDS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    ":",
]
AMINO_ACID_THREE_TO_ONE = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}
PSSM_AA_ORDER = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "X",
]
BLOSUM_DF = pd.read_csv(DATA_DIR / "blosum" / "blosum62.csv").set_index("aa")
BLOSUM_MATRIX = {
    aa_one: {aa_two: BLOSUM_DF.loc[aa_one, aa_two] for aa_two in CANONICAL_AMINO_ACIDS}
    for aa_one in CANONICAL_AMINO_ACIDS
}

MPNN_HIT_TYPE_ORDER = ["score", "recovery", "mutation", "frequency"]
HYDROPHOBIC_ELEMENTS = ["C", "S"]

PLDDT_NORM_FACTOR = 100
PAE_NORM_FACTOR = 10
