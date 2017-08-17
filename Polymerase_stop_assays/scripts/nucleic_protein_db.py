#!/opt/local/bin/python2.7

__doc__ = """
Keep dictionaries of which residues are proteins, which NA, etc...

"""

nuc_resnm = {
    "DA",  "DA3", "DA5", "DAB", "DALA", "DAN", "DC",  "DC3", "DC5", "DCN",
    "DG",  "DG3", "DG5", "DGN", "DT",   "DT3", "DT5", "DTN", "RA",  "RA3",
    "RA5", "RAN", "RC",  "RC3", "RC5",  "RCN", "RG",  "RG3", "RG5", "RGN",
    "RU",  "RU3", "RU5", "RUN", "NA5",  "NA3", "NA",  "NG5", "NG3", "NG",
    "NC5", "NC3", "NC",  "NU5", "NU3",  "NU",  "NT5", "NT3", "NT",  "FA5",
    "FA3", "FA",  "FG5", "FG3", "FG",   "FC5", "FC3", "FC",  "FU5", "FU3",
    "FU",  "FT5", "FT3", "FT",  "A",    "T",   "U",   "C",   "G",   "AMO",
    "GMO", "TMO", "CMO", "PMO", "DP5",  "DP",  "DP3", "URA", "ADE", "CYT",
    "GUA", "THY", "UR2", "AD2"}

protein_dict = {
    'CYS': 'C',
    'ASP': 'D',
    'SER': 'S',
    'GLN': 'Q',
    'LYS': 'K',
    'ILE': 'I',
    'PRO': 'P',
    'THR': 'T',
    'PHE': 'F',
    'ASN': 'N',
    'GLY': 'G',
    'HIS': 'H',
    'LEU': 'L',
    'ARG': 'R',
    'TRP': 'W',
    'ALA': 'A',
    'VAL': 'V',
    'GLU': 'E',
    'TYR': 'Y',
    'MET': 'M'
}
nucleic_dna_dict = {'CYT': 'DC', 'THY': 'DT', 'ADE': 'DA', 'GUA': 'DG'}
nucleic_rna_dict = {'CYT': 'RC', 'URA': 'RU', 'ADE': 'RA', 'GUA': 'RG'}

nucleic_type_dict = {
    'DC': 'py',
    'DC3': 'py',
    'DC5': 'py',
    'DT': 'py',
    'DT3': 'py',
    'DM': 'py',
    'DM3': 'py',
    'DM5': 'py',
    'DH': 'py',
    'DH3': 'py',
    'DH5': 'py',
    'DT5': 'py',
    'DG': 'pu',
    'DG3': 'pu',
    'DG5': 'pu',
    'DA': 'pu',
    'DA3': 'pu',
    'DA5': 'pu',
    'C': 'py',
    'T': 'py',
    'G': 'pu',
    'A': 'pu'
}

py = {'N9': 'N1', 'C4': 'C2', 'N3': 'O2', 'C5': 'N3'}

na_bb_atoms = [
    "P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C1'", "C3'", "C2'", "O3'"
]
