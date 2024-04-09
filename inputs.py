# Input parameters:
# output_dir = f"../results/output/" CodeOcean users.
output_dir = f"output/"
mixtures_needed = 1
density = 10.0
layer_offset = 5.0
number_of_trials = 500
energy_threshold = 0.0
cutoff_distance = 10.0

# Molecules Dictionary:
molecules_dictionary = {
    "benzene": {
        "nbr_of_mols": 500,
        "smiles": "CC(C)Cc1cc2c3c(cc4c(CC(CC4CC)CCCC)c3c1)c1cc(O)cc3c1c2cc(CCC)c3CCC(C)C",
        "rotate": True,
    },
}

# Constants and Unit Conversions:
AVOGADRO_NUMBER = 6.022e23
