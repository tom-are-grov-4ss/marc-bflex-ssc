import logging

import pandas as pd

log = logging.getLogger(__name__)


def mpf_reader(path):

    '''
    Reads the values in 'file.mpf' AS DATAFRAME.

    Input:
        path : Full path to file, including '.mpf' ending

    Returns:
        df : dataframe with values that are read from '.mpf'-file
             {type:dataframe}

    TODO:
    - Consistent in returning either list or np.array
    - Fix df.columns line
    '''

    # Read legends
    with path.open() as f:
        for line in f:
            if 'LEGENDS' in line:
                # Remove ' and split into list, too many spaces and legends
                legends = [leg.replace('\'', '') for leg in line.split('\'\'')]
                legends = [' '.join(leg.split()) for leg in legends]
                legends = [leg.replace('LEGENDS ', '') for leg in legends]

    # Read data into pandas dataframe
    df = pd.read_csv(path, skiprows=8, sep=r'\s+', engine='python')
    df.iloc[0, :-1] = df.iloc[0, 1:].values  # shift first row 1 left
    df = df.dropna(axis=1, how='all')  # drop nas

    legends = legends[:len(df.columns)]

    # Rename columns
    if 'elmom' in path.stem:
        # Bflex bug - not all column names are included
        df.columns = ['Loadstep'] + [f'{c}{n}'
                                     for n in range(1, int(((len(df.columns)
                                                             - 1)/3)) + 1)
                                     for c in ['CORE', 'TEND', 'TOTL']]
    #     elif 'elcurv' in path:
    #         df.columns = [f'Pos_{i/10}' for i in (range(0, 210))]
    else:
        df.columns = legends

    df.drop(df.columns[[0]], axis=1, inplace=True)  # drop column with loadstep

    return df