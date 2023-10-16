#!/usr/bin/env python

from ConfigParser import ConfigParser
from argparse import ArgumentParser
from tracks import trackData


def get_columns(configfile, state=None, cat2bool=True, hcdiff=False):
    conf = ConfigParser()
    conf.read(configfile)
    config = {'Tracks': conf.items('Tracks')}

    ### TRACK PREPARATION ###
    num = 0
    columnNames = []
    columnTypes = []
    if state is not None:
        columnTypes = [bool]
        columnNames = ['y']
        num += 1
    indicatorNames = []
    for track, status in config['Tracks']:
        # check track availability
        assert (track in trackData.keys())

        if status.lower() in ['d', 'default']:
            status = trackData[track]['default']

        # ignore unwanted tracks and remove scaling from unscaled ones
        if status.lower() in ['s', 'scale', 'scaling']:
            pass  # keep scaling on
        elif status.lower() in ['true', 't', 'y', 'yes']:
            if hasattr(trackData[track], 'scaling'):
                del trackData[track]['scaling']
        else:
            continue

        if 'colname' in trackData[track].keys():
            trackData[trackData[track]['colname']] = trackData[track]
            track = trackData[track]['colname']

        # assign track indicies so that combined tracks can assess these
        trackData[track]['id'] = num
        if cat2bool and trackData[track]['type'] is list:
            num += len(trackData[track]['categories'])
        elif trackData[track]['type'] == 'combined':
            trackData[track]['baseLen'] = len(trackData[trackData[track]['base']]['categories'])
            trackData[track]['childLen'] = len(trackData[track]['child'])
            num += trackData[track]['baseLen'] * trackData[track]['childLen']
        else:
            num += 1

        # enforce storing of bools as 0 or 1
        if trackData[track]['type'] is bool:
            trackData[track]['type'] = int

        # generate column names
        if trackData[track]['type'] == 'combined':
            baseNames = trackData[trackData[track]['base']]['categories']
            childNames = trackData[track]['child']
            # child tracks need to be configured before the parent track
            childNames = [child for child in childNames if child in columnNames]
            # remove combinations that are not there (like ignored ones)
            trackData[track]['child'] = childNames
            columnNames.extend([str((base, child)) for child in childNames for base in baseNames])
            columnTypes.extend([trackData[child]['type'] for child in childNames for base in baseNames])
        elif cat2bool and trackData[track]['type'] is list:
            columnNames.extend([str((track, cat)) for cat in trackData[track]['categories']])
            columnTypes.extend([bool for cat in trackData[track]['categories']])
        else:
            columnNames.append(track)
            columnTypes.append(trackData[track]['type'])
        if 'indicator' in trackData[track].keys():
            indicatorNames.append(track + '.ind')

        # do some replacements incase ref and alt are mixed
        if hcdiff:
            for old, rep in [('colname', 'hcdiff_colname'),
                             ('derive', 'hcdiff_derive'),
                             ('transformation', 'hcdiff_transformation')]:
                if rep in trackData[track].keys():
                    trackData[track][old] = trackData[track][rep]

    columnNames.extend(indicatorNames)
    columnTypes.extend([bool for i in indicatorNames])

    return (columnNames, columnTypes, trackData, config)


if __name__ == '__main__':
    parser = ArgumentParser(description="%prog name")
    parser.add_argument("-o", "--output", dest="output", type=str,
                        required=True,
                        help="Name the generated file")
    parser.add_argument("-c", "--config", dest="config", type=str,
                        default='lib/exampleTracks.cfg',
                        help="Config file that specifies used tracks")
    parser.add_argument("-y", dest="state", action='store_true',
                        help="Associated y-Value (default:off)")
    parser.add_argument("-b", "--cat2bool", dest="cat2bool", action='store_true',
                        help="Specify whether categories are split into multiple boolean classifier")

    args = parser.parse_args()

    columnNames, columnTypes, _, _ = get_columns(args.config,
                                                 args.state,
                                                 args.cat2bool)

    columns = zip(columnNames, [str(t) for t in columnTypes])
    with open(args.output, 'w') as stdout:
        stdout.write('\n'.join(['\t'.join(f) for f in columns]))
