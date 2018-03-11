# -*- coding: utf-8 -*-

"""
energy.py - functions for handling energy database
"""

import os
import json
import requests

from collections import defaultdict
from entropy import compute_entropy,compute_boltzmann

# This should remain the same.  Represents amino character mapping to
# common codes used inside FoldX.
codes = {"ALA": "A",
         "ARG": "R",
         "ASN": "N",
         "ASP": "D",
         "ASX": "B",
         "CYS": "C",
         "GLU": "E",
         "GLN": "Q",
         "GLX": "Z",
         "GLY": "G",
         "HIS": "H",
         "ILE": "I",
         "LEU": "L",
         "LYS": "K",
         "MET": "M",
         "PHE": "F",
         "PRO": "P",
         "SER": "S",
         "THR": "T",
         "TRP": "W",
         "TYR": "Y",
         "VAL": "V",
         "TPO": "X"}

def load_db(databaseurl):
    """Either read cache or download the database"""
    localfn = os.path.basename(databaseurl)
    if os.path.exists(databaseurl):
        localfn = databaseurl
    if not os.path.exists(localfn):
        data = requests.get(databaseurl, timeout=30)
        with open(localfn, 'w') as outputfile:
            outputfile.write(data.text)

    with open(localfn, 'rb') as local:
        for line in local:
            yield json.loads(line.rstrip())

def combine_energy_mutations(energydb):
    """This reads line-by-line energy database entries and collapses the
    matching site,wild_type -> mutation mappings to a single
    site,wild_type -> mutations list mapping"""
    energies = defaultdict(dict)
    for jsl in energydb:
        for k,entry in jsl.items():
            protein = entry['protein']
            subprotein = entry['subprotein']
            wt = entry['wt']
            mut = entry['mutation']
            site = entry['site']
            chains = entry['chains']
            epitope = entry.get('epitope', '')
            peptide = entry.get('peptide', '')
            energy = entry['energy_deltas']['total energy']
            # The key needs to be unique enough to match different chains and
            # epitope peptides, but not duplicate entries on top of each other
            key = '{},{},{},{},{},{}'.format(protein, site, wt, epitope,
                                             peptide, chains)
            energies[key]['protein'] = protein
            energies[key]['subprotein'] = subprotein
            energies[key]['epitope'] = epitope
            energies[key]['peptide'] = peptide
            energies[key]['site'] = site
            energies[key]['chains'] = chains
            energies[key]['wt'] = wt
            energies[key]['start'] = entry.get('start', '')
            energies[key]['end'] = entry.get('end', '')
            energies[key][mut] = energy
    return energies

def add_entropies(energydb):
    """For a given energy database calculated in energy.py, calculate
    the entropy for the site"""
    for k in energydb.keys():
        energies = [energydb[k][amino] for amino in codes.values()]
        boltzman_probs = compute_boltzmann(energies)
        # TODO: we could add the probabilities to the dict here...
        energydb[k]['shannon_entropy'] = compute_entropy(boltzman_probs)
    return energydb
