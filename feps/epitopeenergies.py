"""
main.py - main functionality for feps tool
"""
from __future__ import print_function

import click
import csv
import sys

from collections import defaultdict

from feps.energy import load_db,combine_energy_mutations,codes,add_entropies
from feps.entropy import mean

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_epitope_energies(energies):
    """Create a map of epitopes -> epitopeinfo, [wt -> energies]
    Users can then iterate via the peptide chain to grab the energies in order

    This includes the structural energie, which is the mean of the
    site entropies
    """
    epitopedb = defaultdict(dict)
    for _,entry in energies.items():
        # Keys for the epitope db are unique to the epitope, so:
        #   name, peptide, startsite, endsite
        name = entry['epitope']
        peptide = entry['peptide']
        startsite = entry['start']
        endsite = entry['end']
        site = entry['site']
        if not name:
            eprint('Skipping site {} as it does not match an epitope'.format(
                site))
            continue
        key = '{},{},{},{}'.format(name, peptide, startsite, endsite)

        epitopedb[key]['epitope'] = name
        epitopedb[key]['peptide'] = peptide
        epitopedb[key]['peptide_status'] = entry['peptide_status']
        epitopedb[key]['startsite'] = int(startsite)
        epitopedb[key]['endsite'] = int(endsite)
        epitopedb[key]['protein'] = entry['protein']

        wt = entry['wt']
        peptide_state = list(epitopedb[key].get('peptide_state',
                                                ("-" * len(peptide))))
        idx = int(site)-int(startsite)
        if peptide_state[idx] != '-':
            eprint('Overwriting peptide state {} at idx {}!'.format(
                ''.join(peptide_state), idx))
        if peptide[idx] != wt:
            eprint('Peptide mismatch at {}: {} expected, got {}'.format(
                idx, peptide[idx], wt))
        peptide_state[idx] = wt
        epitopedb[key]['peptide_state'] = ''.join(peptide_state)

        # Average energy for the wt at this peptide index
        wtk = '{},{}'.format(idx,wt)
        epitopedb[key][wtk] = mean([entry[x]
                                    for x in codes.values()
                                    if x in entry])
        epitopedb[key][wtk + '-entropy'] = entry['shannon_entropy']

    # Now average energy and structural entropy
    for _,v in epitopedb.items():
        ps = v['peptide_state']
        keys = ['{},{}'.format(i,ps[i]) for i in range(len(ps))
                if ps[i] != '-']
        v['structural_entropy'] = mean([v[k + '-entropy'] for k in keys])
        v['average_energy'] = mean([v[k] for k in keys])

    return epitopedb

@click.command()
@click.option('--database',
              default='https://epitopedata.flowpharma.com/EpitopeData.json',
              help='URL to a jsonl encoded file to dump')
@click.option('--ignore-mutation', default=[], multiple=True,
              type=click.Choice(codes.values()),
              help='Ignore these mutations')
def epitope_energies(database, ignore_mutation):
    db = load_db(database)
    amino_codes = [aa for aa in codes.values() if aa not in ignore_mutation]

    energies = combine_energy_mutations(db, amino_codes)
    energies = add_entropies(energies, amino_codes)

    epitope_energies = get_epitope_energies(energies)

    sorted_keys = sorted(epitope_energies.keys(),
                         key=lambda x: int(x.split(',')[2]))
    fieldnames = [ 'protein', 'epitope', 'peptide', 'peptide_status',
                   'peptide_state',
                   'startsite', 'endsite', 'average_energy',
                   'structural_entropy' ]
    writer = None
    for k in sorted_keys:
        v = epitope_energies[k]

        # Print header
        if not writer:
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames,
                                    extrasaction='ignore')
            writer.writeheader()

        writer.writerow(v)


if __name__ == '__main__':
    epitope_energies()
