"""
main.py - main functionality for tsgen tool
"""
from __future__ import print_function

import click
import csv
import json
import os
import sys
import requests

from collections import defaultdict

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

proteins = [ 'RT', 'TAT', 'P24', 'INT', 'PRO', 'P17', 'REV',
             'GP120', 'NEF' ]

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
    '''Either read cache or download the database'''
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

@click.command()
@click.option('--database',
              default='https://epitopedata.flowpharma.com/EpitopeData.json',
              help='URL to a jsonl encoded file to dump')
def cli(database):
    db = load_db(database)

    energies = defaultdict(dict)
    for jsl in db:
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
            key = '{},{},{},{},{},{}'.format(protein, site, wt, epitope,
                                             peptide, chains)
            energies[key]['protein'] = protein
            energies[key]['subprotein'] = subprotein
            energies[key]['epitope'] = epitope
            energies[key]['peptide'] = peptide
            energies[key]['site'] = site
            energies[key]['chains'] = chains
            energies[key]['wt'] = wt
            energies[key][mut] = energy

    sorted_keys = sorted(energies.keys(), key=lambda x: int(x.split(',')[1]))
    fieldnames = [ 'protein', 'subprotein', 'epitope', 'peptide', 'site',
                   'chains', 'wt' ]
    fieldnames.extend(sorted(codes.values()))
    writer = None
    for k in sorted_keys:
        v = energies[k]

        # Check for missing values and print info to create job to regenerate
        for c in codes.values():
            if c not in v:
                eprint('Missing Data: {}, {}, {}, {}'.format(v['protein'],
                                                             v['site'],
                                                             v['wt'],
                                                             c))

        # Print header
        if not writer:
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
            writer.writeheader()

        writer.writerow(v)


if __name__ == '__main__':
    cli()
