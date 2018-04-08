# -*- coding: utf-8 -*-
"""
main.py - main functionality for feps tool
"""
from __future__ import print_function

import click
import csv
import sys

from feps.energy import load_db,combine_energy_mutations,codes,add_entropies

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

@click.command()
@click.option('--database',
              default='https://epitopedata.flowpharma.com/EpitopeData.json',
              help='URL to a jsonl encoded file to dump')
@click.option('--entropy/--no-entropy', default=False,
              help='Add shannon entropy to outputs')
@click.option('--ignore-mutation', default=[], multiple=True,
              type=click.Choice(codes.values()),
              help='Ignore these mutations')
def cli(database, entropy, ignore_mutation):
    db = load_db(database)
    amino_codes = [aa for aa in codes.values() if aa not in ignore_mutation]

    energies = combine_energy_mutations(db, amino_codes)
    sorted_keys = sorted(energies.keys(), key=lambda x: int(x.split(',')[1]))
    fieldnames = [ 'protein', 'subprotein', 'epitope', 'peptide',
                   'peptide_status', 'site', 'chains', 'wt' ]
    if entropy:
        energies = add_entropies(energies, amino_codes)
        fieldnames.insert(5, 'shannon_entropy')

    fieldnames.extend(sorted(codes.values()))
    writer = None
    for k in sorted_keys:
        v = energies[k]

        # Check for missing values and print info to create job to regenerate
        for c in amino_codes:
            if c not in v:
                eprint('Missing Data: {}, {}, {}, {}'.format(v['protein'],
                                                             v['site'],
                                                             v['wt'],
                                                             c))

        # Print header
        if not writer:
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames,
                                    extrasaction='ignore')
            writer.writeheader()

        writer.writerow(v)


if __name__ == '__main__':
    cli()
