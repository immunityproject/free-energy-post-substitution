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

def get_energy_vectors(db, protein, mutation_protein, start, end):
    start = int(start)
    end = int(end)
    if mutation_protein not in db[protein]:
        raise IOError('Could not find {}, {}, {}'.format(protein, start, end))

    value = db[protein][mutation_protein]
    sites = sorted([int(s) for s in value['energies'].get(protein, {}).keys()])
    if not sites or start not in sites or end not in sites:
        raise IOError('Sites {}, {} are not in {}, {}'.format(start, end,
                                                              protein,
                                                              mutation_protein))

    return {k: v for k, v in value['energies'][protein].items()
            if int(k) >= start and int(k) <= end}

def dump_json(protein, mutation_protein, startsite, endsite, evs):
    print('EVs for protein {}, mutation {}, start {}, end {}: '.format(
        protein, mutation_protein, startsite, endsite))
    print(json.dumps(evs, indent=2))

def dump_csv(protein, mutation_protein, evs):
    fieldnames = [ 'protein', 'mutation_protein', 'site', 'wt' ]
    fieldnames.extend([x['mutation'] for x in list(evs.values())[0]])
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    for site,evlist in evs.items():
        row = dict()
        row['protein'] = protein
        row['mutation_protein'] = mutation_protein
        row['site'] = site
        row['wt'] = evlist[0]['wt']
        for ev in evlist:
            row[ev['mutation']] = ev['energyDelta']
        writer.writerow(row)

def dump_protein(db, protein, mutation_proteins, startsite, endsite, json, csv):
    # If mutation protein is set, check valid and that start/end are not set
    print_sites = []

    if startsite and endsite:
        print_sites.append((startsite, endsite))
    else:
        for mp in mutation_proteins:
            # TODO: Error messages
            if mp not in db[protein].keys():
                print('Mutation protein {} not in protein {}'.format(mp,
                                                                     protein),
                      file=sys.stderr)
                continue

            if not db[protein][mp]['energies']:
                print('Energy lists for protein {} mutation protein {} is'
                      ' empty!'.format(protein, mp), file=sys.stderr)
                continue

            sites = [int(x)
                     for x in db[protein][mp]['energies'][protein].keys()]
            startsite = min(sites)
            endsite = max(sites)

            print('MP: {}, {}, {}'.format(mp, startsite, endsite),
                  file=sys.stderr)

            print_sites.append((mp, str(startsite), str(endsite)))

    for mp,start,end in print_sites:
        evs = get_energy_vectors(db, protein, mp, start, end)
        if not evs or not evs.values():
            print('Could not get enery vectors for protein {}, {},'
                  ' {}, {}'.format(protein, mp, start, end), file=sys.stderr)
        if json:
            dump_json(protein, mp, startsite, endsite, evs)
        if csv:
            dump_csv(protein, mp, evs)

@click.command()
@click.option('--database', default='https://epitopedata.flowpharma.com/P17',
              help='The URL to the epitope data')
@click.option('--protein', '-p', type=click.Choice(proteins), multiple=True,
              help='The selected HIV Protein')
@click.option('--mutation-protein', '-m', default=None, multiple=True,
              help=('The mutation protein, will return the 20 delta G '
                    'energy vectors for all sites'))
@click.option('--json/--no-json', default=False, help='Dump json output')
@click.option('--csv/--no-csv', default=True, help='Dump csv output')
@click.argument('startsite', default=None, type=int, required=False)
@click.argument('endsite', default=None, type=int, required=False)
def cli(database, protein, mutation_protein, startsite, endsite, json, csv):
    db = load_db(database)

    # Select all proteins if none are selected
    if not protein:
        protein = list(db.keys())

    # select all mutation proteins if none are selected
    if not mutation_protein:
        mutation_protein = list()
        for p in protein:
            mutation_protein.extend(list(db[p].keys()))
    # Remove duplicates
    mutation_protein = list(set(mutation_protein))

    for p in protein:
        dump_protein(db, p, mutation_protein, startsite, endsite, json, csv)


if __name__ == '__main__':
    cli()
