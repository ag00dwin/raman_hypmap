def bgr_process(_df,location):

    '''
    MODULE: BACKGROUND REDUCTION (BGR)

    iterative polynomial smoothing method
    As described: Wang, T. and Dai, L. Background Subtraction of Raman
    Spectra Based on Iterative Polynomial Smoothing. Applied Spectroscopy
    2016;71:1169-1179, doi: 10.1177/0003702816670915

    current settings [can be updated within ''filter step'']
    -----
        iteration_accuracy  = 10**-2     # itteration accuracy  
        window_width        = 151               # filter window width 
        polynomial_order    = 2                 # polynomial order 

    prerequisites
    -----
        numpy
        scipy.signal
        sys

    input
    -----
    ''_df'' is dataframe organised with:
        - columns 0 and 1 holding x and y information
        - columns 2: holding the spectral information 
        - row 0 holding the raman shift values for each spectra sample point
    ''location'' is the folder to export the generated csv files

    output
    -----
        the 'bgr' function will automatically save two .csv files to the given location: 
            _bgr.csv; a file in the same format with each spectra background reduced
            _bg.csv; a file in the same format with the background subtracted from each spectra
        the function will also return a pandas dataframe identicle to the _bgr.csv 

    progress bar is thanks to user: https://stackoverflow.com/users/320726/6502
    full details: https://stackoverflow.com/questions/6169217/replace-console-output-in-python

    '''
    
    import numpy as np
    import pandas as pd
    import scipy.signal as sig
    import sys

    # define progress bar functions
    def start_progress(title):
        global progress_x
        sys.stdout.write(title + ": [" + "-"*40 + "]" + chr(8)*41)
        sys.stdout.flush()
        progress_x = 0
        return None
    def progress(x):
        global progress_x
        x = int(x * 40 // 100)
        sys.stdout.write("#" * (x - progress_x))
        sys.stdout.flush()
        progress_x = x
        return None
    def end_progress():
        sys.stdout.write("#" * (40 - progress_x) + "]\n")
        sys.stdout.flush()
        return None
    
    # define filter
    def filter(y_input):          

        global x, y, eratio, iteration_accuracy, sg_background
        # set variables
        iteration_accuracy  = float(10**-2)     # itteration accuracy (IA)  
        window_width        = 151               # filter window width 
        polynomial_order    = 2                 # polynomial order (fixed)

        # apply Savitzky-Golay smoothing filter 
        sg_background = sig.savgol_filter(y_input, window_width, polynomial_order)  
        sg_y          = y_initial-sg_background

        # create comparative index
        eratio = (abs(min(sg_y)))/(abs(max(sg_y)))
        eratio = "{:.6f}".format(eratio)
        eratio = float(eratio)

    # set origin to zero
    _df[0] = _df[0].sub(_df[0].min())
    _df[1] = _df[1].sub(_df[1].min())

    start_progress("Background Reduction")
    total_len = len(_df)-1

    raw_spectra = []
    bgr_spectra = []
    bg_spectra  = []
    # exstract each row of spectra intensity values from the dataframe 
    for index, row in _df.iloc[1:].iterrows():
        percentage_completion = (index/total_len)*100
        progress(percentage_completion)
        
        spectra = np.array(row[2:],dtype='float')
        spectra[np.isnan(spectra)] = 0  # remove nan

        y_initial   = spectra  # set initial pre-BG corrected spectrum
        y_input     = spectra  # set first input spectrum for BG correction
        
        
        # apply filter
        filter(y_input)
        # compare with itteraciton accuracy and apply filter intil accurcy is reached
        while (eratio > iteration_accuracy):          
            filter(np.minimum(y_initial,sg_background)) 
        
        bgr = np.array(y_initial-sg_background,dtype=float)
        bg  = np.array(sg_background,dtype=float)

        # append to array 
        
        raw_spectra.append(np.array(spectra,dtype=float))
        bgr_spectra.append(bgr)
        bg_spectra.append(bg)
        
    end_progress()

    # create new dataframes with background values
    _raw = pd.DataFrame(raw_spectra)    
    _bgr = pd.DataFrame(bgr_spectra)    # background reduced spectra
    _bg  = pd.DataFrame(bg_spectra )    # background that was deducted 
    # duplicate original dataframe and replace spectra with values
    _df_raw = _df.copy()
    _df_bgr = _df.copy()
    _df_bg  = _df.copy()
    _df_raw.iloc[1: , 2:] = _raw
    _df_bgr.iloc[1: , 2:] = _bgr
    _df_bg.iloc[1: , 2:]  = _bg
    # write to seperate csv files and save
    _df_raw.to_csv(location+'_raw.csv',index=False, header=False) # export to csv
    _df_bgr.to_csv(location+'_bgr.csv',index=False, header=False) # export to csv
    _df_bg.to_csv( location+'_bg.csv' ,index=False, header=False) # export to csv
    # create file location

    print('BGR complete')
    print('BGR file saved: '+location+'_bgr.csv')
    print('BG file saved: ' +location+'_bg.csv' )

    #return _df_bgr

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

def hypmap_plot(_df_filter,_df_spectra,file):

    '''
    MODULE: INTERACTIVE MAP 

    creates maps where masks can be drawn to average and save spectra from the supplied dataframe
    maps are created from previously filtered spectra windows

    prerequisites
    -----
        numpy
        pandas
        seaborn

        matplotlib.pyplot
        matplotlib.widgets
        matplotlib.path
        matplotlib.ticker

    input
    -----
    ''_df_filter'' is dataframe organised with:
        - row 0 holding the names for each filter
        - column 0 and 1 holding x and y information
        - column 2-n holding the intensity values for each phase
    ''_df_spectra'' is dataframe organised with:
        - columns 0 and 1 holding x and y information
        - columns 2: holding the spectral information 
        - row 0 holding the raman shift values for each spectra sample point
    ''location'' is the folder to export any saved spectra and related image files

    output
    -----
        the interactive map will show the a plot of the selected spectra points. Clicking ENTER when
        this has loaded will automatically save the X and Y data into a .txt file as well as a image 
        of the map as seen in it's current state. 
            these files are saved as:
                _spectra.txt, organised with:
                    - column 0; r-shift values
                    - column 1; spectral intensity values
                    - row 0; column names 
                _spectra.png, an image of the hyperspectral map at the state when enter is pressed

    '''

    import numpy as np
    import pandas as pd
    import seaborn as sns

    import matplotlib.pyplot as plt
    from matplotlib.widgets import LassoSelector
    from matplotlib.path import Path
    from matplotlib.ticker import FormatStrFormatter

    # select hyperspectral pixels from a matplotlib collection using `LassoSelector`.
    class SelectFromCollection:
        """
        Selected pixels are saved in the `ind` attribute.
        from https://matplotlib.org/stable/gallery/widgets/polygon_selector_demo.html

        Parameters
        -----
        ax : '~matplotlib.axes.Axes'                                ; Axes to interact with.
        collection : 'matplotlib.collections.Collection' subclass   ; collection you want to select from.
        alpha_other : 0 <= float <= 1:
            - To highlight a selection, this tool sets all selected points to an
            - alpha value of 1 and non-selected points to *alpha_other*.
        """

        def __init__(self, ax, collection, alpha_other=0):
            
            self.canvas      = ax.figure.canvas
            self.collection  = collection
            self.alpha_other = alpha_other

            self.xys = collection.get_offsets()
            self.Npts = len(self.xys)

            # ensure that there are separate colors for each object
            self.fc = collection.get_facecolors()
            if len(self.fc) == 0:
                raise ValueError('Collection must have a facecolor')
            elif len(self.fc) == 1:
                self.fc = np.tile(self.fc, (self.Npts, 1))

            self.lasso = LassoSelector(ax, onselect=self.onselect)
            self.ind = []

        def onselect(self, verts):
            path = Path(verts)
            self.ind = np.nonzero(path.contains_points(self.xys))[0]
            self.fc[:, -1] = self.alpha_other
            self.fc[self.ind, -1] = 1
            self.collection.set_facecolors(self.fc)
            self.canvas.draw_idle()

        def disconnect(self):
            self.lasso.disconnect_events()
            self.fc[:, -1] = 1
            self.collection.set_facecolors(self.fc)
            self.canvas.draw_idle()

    # exstract r_shift values from dataframe
    for index, row in _df_spectra.iloc[0:1].iterrows():
        R = np.array(row[2:],dtype='float')
        R[np.isnan(R)] = 0 # remove nan

    # exstract spectral intensity values and append to rows of numpy array
    all_spectra = []
    for index, row in _df_spectra.iloc[2:].iterrows():    
        spectra = np.array(row[2:],dtype='float')
        spectra[np.isnan(spectra)] = 0
        all_spectra.append(spectra)
    all_spectra = np.array(all_spectra)
    # exstract xy data for plotting on spectral map
    y = _df_filter['Y'].values
    x = _df_filter['X'].values
    
    # define xy data for a scatter plot overtop of the heatmap with points that can be selected 
    y_ = (y*((len(set(y))-1)/np.max(y)))+0.5
    x_ = (x*((len(set(x))-1)/np.max(x)))+0.5

    # plot a hyperspectral heatmap for each filter in '_df_filter' dataframe
    for col in _df_filter.columns[3:]:
        # create new dataframe of values from which to construct map
        _plot_df = pd.DataFrame.from_dict(np.array([_df_filter['X'],_df_filter['Y'],_df_filter[col]]).T)
        _plot_df.columns = ['X_value','Y_value','Z_value']
        _plot_df['Z_value'] = pd.to_numeric(_plot_df['Z_value'])
        pivotted= _plot_df.pivot('Y_value','X_value','Z_value')
    
        fig, ax = plt.subplots()
        ax = sns.heatmap(pivotted,
                                vmin    = 0,
                                vmax    = np.max(_plot_df['Z_value']),
                                annot   = False,
                                cmap    = 'Spectral',
                                square  = True, 
                                                                                )
        select = ax.scatter(x_, y_,   
                                marker  = 's',
                                s       = 20,
                                c       = 'black'
                                                                                )
        
        selector = SelectFromCollection(ax, select)

        def accept(event):
            if event.key == "enter":
                # create mask from sum of all selected spectra
                mask = [] 
                for pixel in selector.ind:
                    _s = all_spectra[pixel]
                    mask.append(_s)
                # apply to each pixel selected in mask
                mask     = np.array (mask)
                mask_avg = np.mean  (mask, axis=0)
                mask_std = np.std   (mask, axis=0)
                mask_med = np.median(mask, axis=0)
                mask_sum = np.sum   (mask, axis=0)
                import scipy.signal as sig
                mask_avg_smoothed = sig.savgol_filter(mask_avg,15,2)

                plt.close(2)
                # close and reopen spectra plot for each mask

                _fig, _ax = plt.subplots()

                _ax.plot        (R, mask_avg_smoothed,
                                    label   ='avg [smoothed]',
                                    zorder  =10
                                            )
                _ax.plot        (R, mask_med,
                                    label   ='med',
                                    zorder  =5
                                            )
                _ax.errorbar    (R, mask_avg, mask_std*2,
                                    color   ='black',
                                    ecolor  ='lightgray',
                                    label   ='avg ±2σ',
                                    zorder  =0
                                            )
                plt.legend()

                plt.title('Mask')
                plt.ylabel('intensity')
                plt.xlabel('raman shift (cm-1)')

                # window geometry for spectra window
                mngr = plt.get_current_fig_manager()
                mngr.window.setGeometry(10,30,600,500)  

                # save text file of spectra
                export1 = np.stack([R, mask_avg, mask_std, mask_avg_smoothed, mask_sum], axis=1)
                np.savetxt(file+'_spectra.txt', export1, delimiter="\t", header="R\tS_avg\tS_std\tS_avg_smooth\tS_sum", comments='')
                # save text file of spectra for fitting
                export2 = np.stack([R, mask_avg], axis=1)
                np.savetxt(file+'_spectra_fitting.txt', export2, delimiter="\t", comments='')
                # save image of map
                fig.savefig(file+'_spectra.png', dpi=400)
                fig.savefig(file+'_spectra.pdf', dpi=400)
                plt.show()
    
                return None

        fig.canvas.mpl_connect("key_press_event", accept)

        # formatting hyperspectral map window
        ax.set_title(col)
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(650,30,600,600)  
        majorFormatter = FormatStrFormatter('%0.1f')
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.set_title(str(col))
        plt.gcf().subplots_adjust(bottom=0.15)
        
        plt.show()

    return None    

