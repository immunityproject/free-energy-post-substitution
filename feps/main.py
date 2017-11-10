"""
main.py - main functionality for tsgen tool
"""

import click
import csv
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

def dump_json(protein, startsite, endsite, evs):
    print('EVs for protein {}, start {}, end {}: '.format(protein,
                                                          startsite, endsite))
    print(json.dumps(evs, indent=2))

def dump_csv(protein, evs):
    fieldnames = [ 'protein', 'site', 'wt' ]
    fieldnames.extend([x['mutation'] for x in list(evs.values())[0]])
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    for site,evlist in evs.items():
        row = dict()
        row['protein'] = protein
        row['site'] = site
        row['wt'] = evlist[0]['wt']
        for ev in evlist:
            row[ev['mutation']] = ev['energyDelta']
        writer.writerow(row)

@click.command()
@click.option('--database', default='https://epitopedata.flowpharma.com/P17',
              help='The URL to the epitope data')
@click.option('--protein', '-p', type=click.Choice(proteins),
              help='The selected HIV Protein')
@click.option('--mutation-protein', '-m', default=None,
              help=('The mutation protein, will return the 20 delta G '
                    'energy vectors for all sites'))
@click.option('--json/--no-json', default=False, help='Dump json output')
@click.option('--csv/--no-csv', default=True, help='Dump csv output')
@click.argument('startsite', default=None, type=int, required=False)
@click.argument('endsite', default=None, type=int, required=False)
def cli(database, protein, mutation_protein, startsite, endsite, json, csv):
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
    if json:
        dump_json(protein, startsite, endsite, evs)
    if csv:
        dump_csv(protein, evs)

if __name__ == '__main__':
    cli()
