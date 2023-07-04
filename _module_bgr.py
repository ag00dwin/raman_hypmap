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

def bgr(_df,location):

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