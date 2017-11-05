"""
main.py - main functionality for tsgen tool
"""

import click
import json
import os
import sys
import requests

proteins = [ 'RT', 'TAT', 'P24', 'INT', 'ZIKA_E', 'PRO', 'P17', 'REV',
             'GP120', 'NEF' ]

def load_db(databaseurl):
    '''Either read cache or download the database'''
    localfn = os.path.basename(databaseurl)
    if not os.path.exists(localfn):
        data = requests.get(databaseurl, timeout=30)
        with open(localfn, 'w') as outputfile:
            outputfile.write(data.text)

    return json.load(open(localfn, 'r'))

def get_energy_vectors(db, protein, start, end):
    start = int(start)
    end = int(end)
    for k,v in db[protein].items():
        sites = [int(s) for s in v['energies'].get(protein, {}).keys()]
        if not sites or start not in sites or end not in sites:
            continue

        return {k: v for k, v in v['energies'][protein].items()
                if int(k) >= start and int(k) <= end}

    raise IOError('Could not find {}, {}, {}'.format(protein, start, end))

@click.command()
@click.option('--database', default='https://epitopedata.flowpharma.com/P17',
              help='The URL to the epitope data')
@click.option('--protein', '-p', type=click.Choice(proteins),
              help='The selected HIV Protein')
@click.option('--mutation-protein', '-m', default=None,
              help=('The mutation protein, will return the 20 delta G '
                    'energy vectors for all sites'))
@click.argument('startsite', default=None, type=int, required=False)
@click.argument('endsite', default=None, type=int, required=False)
def cli(database, protein, mutation_protein, startsite, endsite):
    db = load_db(database)

    # If mutation protein is set, check valid and that start/end are not set
    if mutation_protein:
        # TODO: Error messages
        assert mutation_protein in db[protein].keys()
        sites = db[protein][mutation_protein]['energies'][protein].keys()
        if not startsite:
            startsite = min(sites)
        if not endsite:
            endsite = max(sites)
        assert startsite in sites
        assert endsite in sites

    evs = get_energy_vectors(db, protein, startsite, endsite)

    print('EVs for protein {}, start {}, end {}: '.format(protein,
                                                          startsite, endsite))
    print(json.dumps(evs, indent=2))

if __name__ == '__main__':
    cli()
