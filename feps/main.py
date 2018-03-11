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
def cli(database, entropy):
    db = load_db(database)

    energies = combine_energy_mutations(db)
    sorted_keys = sorted(energies.keys(), key=lambda x: int(x.split(',')[1]))
    fieldnames = [ 'protein', 'subprotein', 'epitope', 'peptide', 'site',
                   'chains', 'wt' ]
    if entropy:
        energies = add_entropies(energies)
        fieldnames.insert(5, 'shannon_entropy')

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
