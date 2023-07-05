def filter_window(_df,location,window_location):

    '''
    MODULE: FILTER WINDOW

    compares spectra to r-shift ranges ('windows') to derive maximum intensity values
    to plot as intensities on hypspectral maps

    prerequisites
    -----
        sys
        os
        numpy
        pandas
        seaborn
        matplotlib.pyplot
        scipy.signal

    input
    -----
    ''_df'' is dataframe organised with:
        - columns 0 and 1 holding x and y information
        - columns 2: holding the spectral information 
        - row 0 holding the raman shift values for each spectra sample point
    ''location'' is the folder to export the generated csv file and images

    ALSO a prexisting file in the same directory ending with '_window_ranges.csv'
        - column 0: phase   ; str name of window
        - column 1: char_min; minimum r-shift value, window open
        - column 2: char_max; maximum r-shift value, window closed
        - column 3: comp_min; minimum r-shift value, window open    [to compare spectral intensities]
        - column 4: comp_max; maximum r-shift value, window closed  [to compare spectral intensities]
        - column 5: i_min   ; minimum intensity of maximum spectral intensity value in window
        - column 6: i_max   ; minimum intensity of maximum spectral intensity value in window

    comp_min & comp_max define the open and close values of a second window. The maximum intensity of the spectra
    in the char_min & char_max filter window must be larger than the maximum intensity in the 'comp# comparison 
    window to qualify. This allows for better filtering of more complex spectrum

    output
    -----
        the 'w_filter' function will automatically save one .csv file to the given location:
            _w_filter.csv; a file organised with:
                - row 0 holding the names for each window taken from _df
                - column 0 and 1 holding x and y information
                - column 2-n holding the intensity values for each phase
        the 'w_filter' function will also automatically save a .png and .pdf file of the hypspectral map 
        to the given location.

    progress bar is thanks to user: https://stackoverflow.com/users/320726/6502
    full details: https://stackoverflow.com/questions/6169217/replace-console-output-in-python

    '''
    
    import sys
    import os
    import numpy as np
    import pandas as pd

    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.signal import savgol_filter

    # define progress bar functions
    def start_progress(title):
        global progress_x
        sys.stdout.write(title + ": [" + "-"*40 + "]" + chr(8)*41)
        sys.stdout.flush()
        progress_x = 0
    def progress(x):
        global progress_x
        x = int(x * 40 // 100)
        sys.stdout.write("#" * (x - progress_x))
        sys.stdout.flush()
        progress_x = x
    def end_progress():
        sys.stdout.write("#" * (40 - progress_x) + "]\n")
        sys.stdout.flush()

    # define compare function
    def mineralcomp(r_shift,intensity,window_df):

        comparison_array = []
        # import windows
        for index, row in window_df.iterrows():
            r_shift = np.array(r_shift)
            # retrieve data for windows from dataframe
            phase    = row['phase']
            char_min = row['char_min']
            char_max = row['char_max']
            comp_min = row['comp_min']
            comp_max = row['comp_max']
            i_min    = row['i_min']
            i_max    = row['i_max']

            # if comparison region is given, collect region of spectra
            r_comp = []
            i_comp = []
            max_i_comp = 0
            if comp_max > 0 :
                for nc in r_shift:                
                    if  comp_min <= nc <= comp_max:
                        r_comp.append(nc)
                for pc in r_comp:
                    i_comp.append(intensity[np.where(r_shift == pc)])
                i_comp     = np.array(i_comp)
                max_i_comp = np.amax(i_comp)

            # collect region of spectra within match range
            r_match = []
            i_match = []
            for ns in r_shift:                
                if  char_min <= ns <= char_max:
                    r_match.append(ns)
            for ps in r_match:
                i_match.append(intensity[np.where(r_shift == ps)])
            i_match     = np.array(i_match)
            max_i_match = np.amax (i_match)    

            ints = 0

            # compare to comparison region 
            # if smaller, set intensity to 0
            if comp_max > 0:
                if max_i_match <= max_i_comp:
                    ints = 0
                else:
                    ints = max_i_match
            else:
                ints = max_i_match

            # ensure within intensity threshold
            # if outside, set intensity to 0
            if max_i_match < i_min:
                ints = 0
            if max_i_match > i_max:
                ints = 0

            # append to array     
            comparison_array.append([phase,ints])

        return comparison_array

    # import dataframe of windows
    path = window_location
    window_df = pd.read_csv(path, sep=',', index_col=False)

    # create empty dataframe to export window comparison data into
    _df_wind_comp = pd.DataFrame()   
    # append x and y columns from input dataframe
    xy = _df.loc[2:,0:1]
    _df_wind_comp = xy.copy()
    _df_wind_comp = _df_wind_comp.rename(columns={_df_wind_comp.columns[0]: 'X'})
    _df_wind_comp = _df_wind_comp.rename(columns={_df_wind_comp.columns[1]: 'Y'})

    # exstract r_shift values from dataframe
    for index, row in _df.iloc[0:1].iterrows():
        R = np.array(row[2:],dtype='float')
        R[np.isnan(R)] = 0 # remove nan

    # exstract spectral intensity values and append to rows of numpy array
    all_spectra = []
    for index, row in _df.iloc[2:].iterrows():    
        spectra = np.array(row[2:],dtype='float')
        spectra[np.isnan(spectra)] = 0
        all_spectra.append(spectra)
    all_spectra = np.array(all_spectra)
    
    # iterate over all spectra and compare with standards
    # create new array to transform into dataframe that contains values
    wind_comp = []
    # add mineral names to first row of array 
    wind_comp.append(np.rot90(mineralcomp(R,all_spectra[0],window_df), k=1, axes=(1, 0))[0])  
    
    # itterate over all spectra and append new row of values for each spectra
    # each spectra behaves as an XY coorindate on hyperspectral map 
    start_progress("Comparing to windows")
    s_len = len(all_spectra)
    for index, s in enumerate(all_spectra):
        percentage_completion = (index/s_len)*100
        progress(percentage_completion)

        output = np.rot90(mineralcomp(R,s,window_df), k=1, axes=(1, 0))
        wind_comp.append(output[1])
    end_progress()

    # reshape columns into individual lists
    for index, item in enumerate(wind_comp[0]):
        name = item
        val = []
        for i, value in enumerate(wind_comp[1:]): 
           val.append(value[index]) 
        val = np.array(val,dtype='float')
        _df_wind_comp[name] = val
    # export to csv
    _df_wind_comp.to_csv(location+'_w_filter.csv')

    # create a hyperspectral image in greyscale 
    # export and save image
    for col in _df_wind_comp.columns[2:]:

        _plot_df = pd.DataFrame.from_dict(np.array([_df_wind_comp['X'],_df_wind_comp['Y'],_df_wind_comp[col]]).T)
        _plot_df.columns = ['X_value','Y_value','Z_value']
        _plot_df['Z_value'] = pd.to_numeric(_plot_df['Z_value'])
        pivotted= _plot_df.pivot('Y_value','X_value','Z_value')

        ax = sns.heatmap(pivotted,
                                vmin    = 0,
                                vmax    = np.max(_plot_df['Z_value']),
                                annot   = False,
                                cmap    = 'Greys_r',
                                square  = True, )
        ax.set_title(str(col))
        # format labels
        from matplotlib.ticker import FormatStrFormatter
        majorFormatter = FormatStrFormatter('%0.1f')
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_major_formatter(majorFormatter)
    
        plt.gcf().subplots_adjust(bottom=0.15)

        figure = ax.get_figure()  
        # export as png and as pdf
        figure.savefig(location+'_png_'+str(col)+'.png', dpi=400)
        figure.savefig(location+'_pdf_'+str(col)+'.pdf', dpi=400)
        # clear plot for next image
        plt.clf()

    # close all plotting windows
    plt.close('all')

    return None


    

