
contrasts = {
        'D-v-C':      ('D-v-C.ALL.D',   'D-v-C.ALL.C'),  # D -v- C 
        'H-v-C':      ('H-v-C.ALL.H',   'H-v-C.ALL.C'),  # H -v- C
        'HD-v-C':     ('HD-v-C.ALL.HD', 'HD-v-C.ALL.C'), # HD -v- C
        'HD-v-H':     ('HD-v-H.ALL.HD', 'HD-v-H.ALL.H'), # HD -v- H
        'W-v-C':      ('W-v-C.ALL.W',   'W-v-C.ALL.C'),  # W -v- C   
}

file_prefix_per_contrast = {
        'D-v-C':    "C-D/network-C_D",
        'H-v-C':    "C-H/network-C_H",
        'HD-v-C':   "C-HD/network-C_HD",
        'HD-v-H':   "H-HD/network-H_HD",
        'W-v-C':    "C-W/network-C_W"
}

file_suffix = "proteomics_metabolomics_hormonomics_qPCR_phenomics.txt"
    
files = {}
for contrast in contrasts:
    files[contrast] = f"{file_prefix_per_contrast[contrast]}-{file_suffix}"