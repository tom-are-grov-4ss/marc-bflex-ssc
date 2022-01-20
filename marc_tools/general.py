import dask
import time
import logging

import numpy as np
import multiprocessing as mp

from dask.distributed import LocalCluster


logging.basicConfig(filename=r'./Marc_tools.log', filemode='w',
                    level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
logging.info('Started')


def dask_client(cpu=0.7, memory_limit="4GB", port='8786'):
    '''
    Function to initiate the dask client.

    Input:
        - cpu : Controlling parallel cores (0.5 is 50 % of all machine cores)
        - memory limit : Memory limit
        - port : port number to use

    Return:
        - client : Dask multiprocessing client
    '''

    client = LocalCluster(
        n_workers=int(cpu * mp.cpu_count()),
        processes=True,
        threads_per_worker=1,
        memory_limit=memory_limit,
        ip=f'tcp://localhost:{port}')

    return client


def readable_time(seconds):
    '''
    Creates a readable string with time info based on seconds

    Input
    - seconds : Seconds [s]

    Return
    - (readable string with time)
    '''
    times = [('year', 365 * 24 * 60 * 60),
             ('day', 24 * 60 * 60),
             ('hour', 60 * 60),
             ('minute', 60),
             ('second', 1)]

    chunks = []
    for name, secs in times:
        qty = seconds // secs
        if qty:
            if qty > 1:
                name += "s"
            chunks.append(f'{int(qty)} {name}')

        seconds = seconds % secs

    try:
        return ', '.join(chunks[:-1]) + ' and ' + chunks[-1] if len(chunks) > 1 else chunks[0]
    except IndexError:
        return '0.0 seconds'


def pretty_round(n):
    '''
    Rounding a (decimal) number to a pretty number

    Input:
        - n : ugly number to be rounded

    Return:
        - n : rounded pretty number
    '''

    if n < 1:
        base = 5
        tens = abs(np.floor(np.log10(n))) + 1
        n = base * np.ceil(n * 10**tens / base) / 10**tens
    else:
        base = 10
        tens = abs(np.floor(np.log10(n)))
        n = 1 / base * np.ceil(n / 10**tens * base) * 10**tens

    return n


def completed(message, tic):
    '''
    Function that writes completed message to log file

    Input:
        - message : message to be displayed
        - tic : start time of process
    '''

    logging.info(f'Completed {message.strip()} in '
                 f'{readable_time(time.time() - tic)}')


def find_nearest(mylist, value):
    '''
    Fix so that it find highest number as well
    '''

    # Sort list
    sorted_list = sorted(mylist)

    if (value < sorted_list[0]) or (value > sorted_list[-1]):
        return False

    # Loop through
    for i, s in enumerate(sorted_list):
        if s <= value:
            pass
        else:
            return [sorted_list[i-1], s]