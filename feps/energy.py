# -*- coding: utf-8 -*-

"""
energy.py - functions for handling energy database
"""

import os
import json
import logging
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

def combine_energy_mutations(energydb, amino_mutations):
    """This reads line-by-line energy database entries and collapses the
    matching site,wild_type -> mutation mappings to a single
    site,wild_type -> mutations list mapping
    """
    energies = defaultdict(dict)
    for jsl in energydb:
        for k,entry in jsl.items():
            mut = entry['mutation']
            if mut not in amino_mutations:
                logging.info('Skipping {} because it is not in the mutations'
                             ' list'.format(k))
                continue
            protein = entry['protein']
            subprotein = entry['subprotein']
            wt = entry['wt']
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
            energies[key]['peptide_status'] = entry.get('peptide_status', '')
            energies[key]['site'] = site
            energies[key]['chains'] = chains
            energies[key]['wt'] = wt
            energies[key]['start'] = entry.get('start', '')
            energies[key]['end'] = entry.get('end', '')
            energies[key][mut] = energy
    return energies

def add_entropies(energydb, aminos, include_absolute_entropy = False,
                  include_boltzman_entropy = False,
                  include_absolute_boltzman_entropy = False):
    """For a given energy database calculated in energy.py, calculate
    the entropy for the site"""
    for k in energydb.keys():
        energies = list()
        for amino in aminos:
            if amino not in energydb[k]:
                logging.error('Missing amino energy from {}: {}'.format(
                    k, amino))
                continue
            energies.append(energydb[k][amino])
        boltzman_probs = compute_boltzmann(energies)
        # TODO: we could add the probabilities to the dict here...
        energydb[k]['shannon_entropy'] = compute_entropy(boltzman_probs)
        if include_absolute_entropy:
            abs_boltzman_probs = compute_boltzmann([abs(e) for e in energies])
            energydb[k]['absolute_shannon_entropy'] = compute_entropy(
                abs_boltzman_probs)
        if include_boltzman_entropy:
            boltzman_classic_probs = compute_boltzmann(
                [(e*4184/((1.38e-23*6.02e+23)*310.15)) for e in energies])
            energydb[k]['boltzman_shannon_entropy'] = compute_entropy(
                boltzman_classic_probs)
        if include_absolute_boltzman_entropy:           
            abs_boltzman_classic_probs = compute_boltzmann(
                [abs(e*4184/((1.38e-23*6.02e+23)*310.15)) for e in energies])
            energydb[k]['absolute_boltzman_shannon_entropy'] = compute_entropy(
                abs_boltzman_classic_probs)
    return energydb
