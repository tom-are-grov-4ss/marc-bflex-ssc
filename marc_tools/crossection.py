# Modules
# import pathlib
import plotly
import logging
import re
import os
import sys

# Modules with short names
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import marc_tools as mt

# Packes from modules
# from IPython.display import display, clear_output, HTML
plotly.offline.init_notebook_mode(connected=True)

log = logging.getLogger(__name__)

# Standard defintions
# TODO : can be written to {'Fx' : {'Force X' : 'Force [N]'}} and refer to ass fm['Fx'].key
fm = {'EI': 'Stiffness', 'Fx': 'Axial Force X', 'Fy': 'Transverse Force Y',
      'Fz': 'Shear Force Z', 'Ftot': 'Total force', 'Mx': 'Torque', 'My': 'Moment Y',
      'Mz': 'Moment Z', 'K1_2': 'Curvature', 'K3_4': 'Curvature',
      'Ovality Y': 'Ovality', 'Iyy': 'Iyy', 'Iyz': 'Iyz', 'Izz': 'Izz',
      'Av Pl Strain': 'Average plastic strain', 'Stress': 'Stress',
      'Max Pl Strain': 'Max plastic strain'}

fm_units = {'EI': 'Stiffness [Nmm2]',  'Fx': 'Force [N]', 'Fy': 'Force [N]',
            'Fz': 'Force [N]', 'Ftot': 'Total Force [N]', 'Mx': 'Moment [Nmm]',
            'My': 'Moment [Nmm]', 'Mz': 'Moment [Nmm]',
            'K1_2': 'Curvature [1/m]', 'K3_4': 'Curvature [1/m]',
            'Ovality Y': 'Ovality [-]', 'Av Pl Strain': 'Plastic strain [-]',
            'Stress': 'Stress (MPa)', 'Iyy': 'Iyy', 'Iyz': 'Iyz', 'Izz': 'Izz',
            'Max Pl Strain': 'Plastic strain [-]'}

# Names used in plot instead of names used in model, e.g. {'model_name' : 'Plot name'}
set_explanation = {'fbr': 'Front brick tooth', 'flh': 'Front liner holder'}

colors_4s = ['#00a0b0', '#f3776f', '#012b5d', '#009ddc',
             '#00a88f', '#faa634', '#c0c1ba', '#f8ed9b',
             '#00a0b0', '#f3776f', '#012b5d', '#009ddc',
             '#00a88f', '#faa634', '#c0c1ba', '#f8ed9b']

line_styles = ['lines', 'lines+markers', 'markers',
               'lines', 'lines+markers', 'markers',
               'lines', 'lines+markers', 'markers']


def cs_reader(file_dir, extra=False, combine=''):
    """
    Reads .csv-file produced by CrossSection
    Renames columns to more readable names
    Calculates "total force"

    Input:
    - file_dir : directory of file (pathLib)
    - combine : combine all bodies matching the names in combine (list)

    TODO: Return von mieses
    """

    # Definitions
    name, new_name = '', ''

    # Specify whether to look after element set or contact body name
    if os.path.isfile(file_dir):
        if 'sum' in file_dir.name:
            cs_type = 'Element Set'
        elif 'cb' in file_dir.name:
            cs_type = 'Contact Body'
        else:
            raise Exception(f'Check {file_dir} name - should contain cb or sum')
    else:
        raise Exception(f'cs_reader did not find {file_dir}')

    # Read csv file
    df = pd.read_csv(file_dir)
    df.dropna(how='all', axis=1, inplace=True)  # removes columns with nan only

    # Combine columns
    df_combined = pd.DataFrame()
    if isinstance(combine, list):
        for layer in combine:
            # Find columns to combine
            df_comb = df[[col for col in df.columns if layer in col]]

            # Find properties
            props = set(col.split('.')[0] for col in df_comb.columns)

            # Iterate over properties and combine columns with match
            for prop in props:
                comb_match = f'{prop}.{layer}'
                df_combined[comb_match] = df_comb[[col for col in
                                                   df_comb.columns if comb_match
                                                   in col]].sum(axis=1)

        # Find columns that have been combined and that are excessive
        columns_final = df_combined.columns
        columns_all = [col for col in df.columns if any([layer in col for layer in combine])]
        columns_excessive = list([col for col in columns_all if col not in columns_final])

        # Remove excessive columns and change content in matched columns
        df.drop(columns_excessive, inplace=True, axis=1)

        # Changes columns
        for col in columns_final:
            df[col] = df_combined[col]

    # Find element data: number, names and order
    for c in df.columns.values:
        if cs_type in c:
            name = df[c][1].strip()
        if name:
            col_name = c.split('.')[0]
            new_name = '.'.join([col_name, name])
            df.rename(columns={c: new_name}, inplace=True)

    # Calc additional parameters for each part
    cbodies = set([n.split('.')[-1] for n in df.columns.values if '.' in n])

    if extra is True:
        for p in sorted(cbodies):
            # Total force
            df['.'.join(['Ftot', p])] = (df['.'.join(['Fx', p])]**2 +
                                        df['.'.join(['Fy', p])]**2 +
                                        df['.'.join(['Fz', p])]**2).pow(1/2)

            # Calc stiffness for each part
            df['.'.join(['EI', p])] = (df['.'.join(['My', p])]
                                    / df['.'.join(['K1_2'])]).replace([np.inf],
                                                                        0.0)

    return df, cbodies


def find_files(run_dir, *args):
    '''
    Function that finds all files with specified ending

    :param run_dir: Directory containing files
    :param args: Arguments which run_dir must contain

    :return files: Files containing arguments

    '''
    files = []

    for arg in args:
        file_paths = sorted(run_dir.rglob(f'*{arg}*'))
        files += [p.relative_to(run_dir) for p in file_paths]

    return files


def find_ranges(my_dict, x_divisor=500, y_divisor=50000):
    '''
    Function that finds the range of the x- and y-axis based on data in my_dict

    :dict my_dict: Dictionary to process
    :param x_divisor: Number the x-range are rounded up to
    :param y_divisor: Number the y-range are rounded up to

    :return x_range, y_range: Two lists with min/max
    '''

    x_max = max(list(my_dict.values())[0].columns.values)
    x_range = [0, ((x_max // x_divisor + 1) * x_divisor)]

    y_max = 0
    y_min = 0
    for layer in my_dict.keys():
        if 'SUM' not in layer:
            for col in my_dict[layer].columns.values:
                my_list = list(my_dict[layer][col])
                my_min = min(my_list)
                my_max = max(my_list)

                if my_min < y_min:
                    y_min = my_min

                if my_max > y_max:
                    y_max = my_max

    y_range_max = ((y_max // y_divisor + 1) * y_divisor)
    y_range_min = ((abs(y_min) // y_divisor + 1) * y_divisor) * -1

    y_range = [y_range_min, y_range_max]

    return x_range, y_range


def read_results(case_path, res):
    '''
    Function that reads all .csv-files in case_path

    :param case_pth: path with .csv-files
    :param res: result to be extracted from csv-file

    TODO:
    - Hva hvis det er .csv-filer som ikke inneholder resultater i mappen?
    '''

    # Initialize vectors
    res_dict = {}
    temp_dict = {}
    # time = pd.DataFrame()
    displacements = pd.DataFrame()
    # temp_df = pd.DataFrame()

    # Find .csv-files with results
    csv_files = find_files(case_path, 'csv')

    # Find names of contact bodies
    df = pd.read_csv(case_path / csv_files[0])  # NEW
    names = [df[cb][0].lstrip() for cb in df.columns.values if 'Contact Body' in cb]

    # Initialize dictionary with dataframe for each layer
    for name in names:
        res_dict[name] = pd.DataFrame()

    # Extract desired result from each file
    for csv_file in csv_files:
        csv_name = csv_file.stem
        node = csv_name.split('_')[-1].replace('cb', '')

        csv_path = case_path / csv_file
        df = pd.read_csv(csv_path)

        temp_dict[node] = extract_results(df, res)

    # Re-format data
    for node in temp_dict.keys():
        displacements[node] = temp_dict[node]['X-disp']

        for name in names:
            res_dict[name][displacements[node][0]] = temp_dict[node][name]

    res_dict['X-disp'] = displacements

    # Sort columns (with forces for axial positions along pipe) numerically for plotting
    for name in res_dict.keys():
        res_dict[name] = res_dict[name].sort_index(axis=1)

    return res_dict


def extract_results(df, res):
    '''
    Function that extract result result from dataframe

    :param df: DataFrame wtih all results
    :param res: Result to be extracted

    :return df_res: DataFrame with selected results
    '''

    if not isinstance(res, str):
        sys.exit('{} is not a string'.format(res))

    df_res = pd.DataFrame()
    # node_disp = pd.DataFrame()

    df_res['Time'] = df['time']
    df_res['X-disp'] = df['Origo X GCS']

    names = [df[cb][0].lstrip() for cb in df.columns.values
             if 'Contact Body' in cb]  # extract names of contact bodies
    res_cols = [col for col in df.columns if res in col]

    if len(res_cols) == len(names):
        for i, col in enumerate(res_cols):
            df_res[names[i]] = df[col]

        return df_res
    else:
        sys.exit('Can not extract {} from dataframe'.format(res))


def trim_list(my_list, rem_strings):
    '''
    Function that removes strings from a list

    :list my_list: List with strings
    :list rem_strings: List with strings to be removed.

    :return my_list: List with strings removed
    '''

    for s in rem_strings:
        if s in my_list:
            my_list.remove(s)
    return my_list


def strain2stress(df):
    '''
    Function used to show stress of material instead of strain,
    which is given by CS as default

    df :

    '''

    project_folder = 'undefined'
    df_res = {}
    # read excel
    mat_316 = pd.ExcelFile(os.path.join(project_folder,
                                        '1324_Material.xlsx')).parse('316')

    for case in df.keys():
        df_stress = pd.DataFrame()
        for col in df[case].columns.values:
            df_stress[col] = np.interp(df[case][col],
                                       mat_316['True plastic strain'],
                                       mat_316['True stress'])

        df_res[case] = df_stress.where(df_stress > mat_316['True stress'][0])

    return df_res


def plot_data(x, y, res, keep='', remove='', plot_title=''):
    """
    Function used to hide default input from main cell

    x : dataframe to be plotted along x-axis
    y : dataframe to be plotted along y-axis
    res : result to be plotted (e.g. 'Fx')
    keep : if to keep only a selected number of sets/cbs
           (e.g. 'brick1', 'brick2')
    remove : if to remove selected sets/cbs (e.g. 'none', 'liner')
    plot_title : title of plot. default is case name
    """

    if not isinstance(res, str):
        raise Exception('Error in plot_data - res not string')

    df_plot = df_filter(y, res, keep, remove)

    data = []
    for i, plot_list in enumerate(df_plot.columns.values):
        name = plot_list.split('.')[-1]
        if name[0:-1] in set_explanation.keys():
            name = set_explanation[name[0:-1]] + name[-1]
        else:
            name = name.capitalize()

        trace = go.Scatter(mode=line_styles[0],
                           x=x,
                           y=df_plot[plot_list],
                           name=name,
                           line={'color': colors_4s[i]}
                           )
        data.append(trace)

    try:
        yaxis_title = fm_units[plot_list.split('.')[0]]
    except (KeyError, UnboundLocalError):
        yaxis_title = 'Plastic strain [-]'
        # yaxis_title = keep

    layout = go.Layout(title=plot_title,
                       xaxis=dict(title='Contact force [N]'),
                       # yaxis=dict(title='Plastic strain [-]'))
                       yaxis=dict(title=yaxis_title))

    plotly.offline.iplot(go.Figure(data=data, layout=layout))


def read_excel_results(excel_file_dir, res, cases):
    '''
    Reads each sheet in excel_file_dir, and only columns matching 'res' will be
    extracted. Sheet name will be changed to readable format based on dictionary
    'cases'

    TODO:
    - What if excel sheet gives multiple matches?
      E.g. if both Contact Force X & Y is present?
    '''
    df_excel = {}
    xl = pd.ExcelFile(excel_file_dir)

    # For each sheet - this gives result on format df_excel[case][result]
    for sheet in xl.sheet_names:
        if any([sheet in case for case in cases.values()]):
            df = xl.parse(sheet)
            df.dropna(how='all', axis=1, inplace=True)

            df_excel[sheet] = df

    return df_excel


def df_filter(df, res='', keep='', remove=''):
    """
    For the dataframe 'df', ONLY the result 'res' will be kept,
    and columns containing words in 'remove' will be removed

    df : dataframe to filter
    res : result to be kept
    keep : strings to be kept
    remove : strings to be removed

    TODO: df.drop lines can def be improved wrt readability
    """

    logging.info(f'{keep}')
    # Check type - converts to lists
    if isinstance(res, str):
        res = [res]

    if isinstance(remove, str):
        remove = [remove]

    if isinstance(keep, str):
        keep = [keep]

    # Operations
    for r in res:
        if r != '':
            df = df.filter(regex=r)

    if keep != ['']:
        # Default keep
        keep = ['inc', 'time', 'K1_2', 'Ovality Y'] + keep

        # Keep
        logging.info(f'{list(df.columns.values)}')
        df = df.drop([i for i in df.columns.values
                      if not any(y for y in keep
                                 if i.lower().endswith(y.lower()))], axis=1)
        logging.info(f'{keep}')

    if remove != ['']:
        # Remove
        df = df.drop([i for i in df.columns.values
                      if any(y for y in remove
                             if i.lower().endswith(y.lower()))], axis=1)

    return df


def df_combine_csv(case, excel_file, post_folder):
    """
    ! U N F I N I S H E D   A N D   U N U S E D !

    For cases where multiple independent cuts and sets have been run,
    this function collects and stores all sets in a dataframe

    case :
    excel_file :
    post_folder :
    """

    df_set = pd.DataFrame()
    # xl = pd.ExcelFile(os.path.join(project_folder, excel_file))

    # # extract info from excel file
    # for sheet in xl.sheet_names:
    #     case_folder = os.path.join(sheet, post_folder)
    #     df_excel = xl.parse(sheet)

    #     if case in sheet:
    #         # read each csv-file and combine into one df
    #         for elset, csv in zip(df_excel['set'], df_excel['csv-file']):
    #             df = df_reader(case_folder, csv)
    #             df_set = pd.concat([df_set, df.filter(regex=elset)], axis=1, sort=False)
    #     else:
    #         print('{} not in {}'.format(case, sheet))

    return df_set


def df_cs(project_folder, runs, exclude_cases=[], positions='start'):
    '''
    Input:
        - project_folder : folder with project {type: pathlib.Path}
        - runs : list with folders with runs
        - exclude_cases : case names to exclude

    Return
        - df_case[case][position][result]
    '''

    # Initialization
    sorted_csv, df_case, increments, cbodies = {}, {}, {}, {}

    # Find csv files
    csv_files = []
    for run in runs:

        run_folder = project_folder / run

        # Find csv files
        all_csv_files = [run_folder / f.lower() for f in os.listdir(run_folder)
                         if len(re.findall(r'\d+cb.csv', f)) == 1]
        csv_files += [f for f in all_csv_files
                      if os.path.isfile(run_folder /
                                        str(f).replace('cb.csv', 'cb.res'))]

    # Remove cases to be exluded
    csv_files = [csv for csv in csv_files
                 if not any([ex.lower() in csv.stem.lower()
                             for ex in exclude_cases])]

    # Sort csv so that csv-files are tied to correct model
    for csv_file in csv_files:
        csv_info = csv_file.stem.rpartition('_')
        model, csv = csv_file.parent / csv_info[0], csv_info[2]

        if model not in sorted_csv.keys():
            sorted_csv[model] = [csv]
        else:
            sorted_csv[model].append(csv)

    # Get dataframe for each case
    for case in sorted_csv.keys():

        # For each crossection file
        df_csv = {}
        for csv_file in sorted_csv[case]:

            # Find directory to csv and read csv file
            csv_file_dir = case.parent / f'{case.stem}_{csv_file}.csv'

            try:
                df_temp, cbodies[case.stem.lower()] = mt.crossection\
                                                        .cs_reader(csv_file_dir)
            except Exception as e:
                log.warning(f'{e} : {csv_file_dir} is empty')
                continue

            # One df for each position
            pos = df_temp['Origo Z GCS'][0]

            # Overwrite coordinates for each position
            if positions == 'start':
                df_temp['Origo Z GCS'] = [df_temp['Origo Z GCS'][0]] * len(df_temp)

            # Put each df into dict
            df_csv[f'{pos:.0f}'] = df_temp
            df_case[case.stem] = df_csv

            # Extract increments
            if case.stem not in increments:
                increments[case.stem] = df_temp['inc'].iloc[-1]

    return df_case
