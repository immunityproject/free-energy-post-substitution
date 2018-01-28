"""
main.py - main functionality for tsgen tool
"""

import click
import csv
import json
import os
import sys
import requests

from collections import defaultdict

proteins = [ 'RT', 'TAT', 'P24', 'INT', 'PRO', 'P17', 'REV',
             'GP120', 'NEF' ]

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
            energy = entry['energy_deltas']['total energy']
            key = '{},{}'.format(protein,site)
            energies[key]['protein'] = protein
            energies[key]['subprotein'] = subprotein
            energies[key]['epitope'] = entry.get('epitope', '')
            energies[key]['peptide'] = entry.get('peptide', '')
            energies[key]['site'] = site
            energies[key]['wt'] = wt
            energies[key][mut] = energy

    sorted_keys = sorted(energies.keys(), key=lambda x: int(x.split(',')[1]))
    fieldnames = [ 'protein', 'subprotein', 'epitope', 'peptide', 'site', 'wt' ]
    cur_fns = set()
    writer = None
    for k in sorted_keys:
        v = energies[k]
        new_fns = {f for f in v.keys()
                   if f not in cur_fns and f not in fieldnames}
        if new_fns or not writer:
            cur_fns = new_fns
            fns = [f for f in fieldnames]
            fns.extend(sorted(new_fns))
            writer = csv.DictWriter(sys.stdout, fieldnames=fns)
            writer.writeheader()
        writer.writerow(v)


if __name__ == '__main__':
    cli()
