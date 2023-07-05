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
                export = np.stack([R, mask_avg, mask_std, mask_avg_smoothed, mask_sum], axis=1)
                np.savetxt(file+'_spectra.txt', export, delimiter="\t", header="R\tS_avg\tS_std\tS_avg_smooth\tS_sum", comments='')
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

