# raman_hypmap

Raman hyperspectral processing via python scripts tailoured for output of Horiba instruments, specifcially the [LabSpec 6 Spectroscopy Suite Software](https://www.horiba.com/int/scientific/products/detail/action/show/Product/labspec-6-spectroscopy-suite-software-1843/) product. Theoretically this code will work on any hyperspectral Raman data if the .txt input is formatted correctly. 

Three python modules are available:
* bgr_process — automated background reduction of entire map
* filter_window — using defined peaks to create spectral intensity maps 
* hypmap_plot — interactive plot of intensity maps to create and save spectral masks

## Referencing

For use of this code, please cite: *In review.* Doi and further information will be added upon publication 

Please also reference Wang and Dai (2016) for the iterative polynomial smoothing method used for background reduction of spectra `bgr_process.py`.

Wang, T. and Dai, L., 2017. Background subtraction of Raman spectra based on iterative polynomial smoothing. Applied spectroscopy, 71(6), pp.1169-1179. doi: https://doi.org/10.1177/0003702816670915

>This algorithm is based on an iterative polynomial smoothing method that highly reduces the need for experience and a priori knowledge. First, a polynomial filter is applied to smooth the input spectrum (the input spectrum is just an original spectrum at the first iteration). The output curve of the filter divides the original spectrum into two parts, top and bottom. Second, a proportion is calculated between the lowest point of the signal in the bottom part and the highest point of the signal in the top part. The proportion is a key index that decides whether to go into a new iteration. If a new iteration is needed, the minimum value between the output curve and the original spectrum forms a new curve that goes into the same filter in the first step and continues as another iteration until no more iteration is needed to finally get the background of the original spectrum. Results from the simulation experiments not only show that the iterative polynomial smoothing algorithm achieves good performance, processing time, cost, and accuracy of recovery, but also prove that the algorithm adapts to different background types and a large signal-to-noise ratio range.

The progress bar used in this code is from [user](https://stackoverflow.com/users/320726/6502) and full details are on [stackoverflow](https://stackoverflow.com/questions/6169217/replace-console-output-in-python)

The pixel selector used for the hyperspectral map `hypmap_plot` is adapted from a [matplotlib tutorial](https://matplotlib.org/stable/gallery/widgets/polygon_selector_demo.html) for creating a polygon selector. 

### License

The project is licensed under the GNU-v3 license.

## Installation

In order to run the package please run the following commands: ???

Alternatively, the files can be downloaded and imported as modules locally. 

## Usage

Import a .txt hyperspectral map file from the horiba software (*e.g.,* [LabSpec 6](https://www.horiba.com/int/scientific/products/detail/action/show/Product/labspec-6-spectroscopy-suite-software-1843/)). The file is structured as a tab deliminated file with the first row denoting the `Raman Shift` (cm<sup>-1</sup>) values. Note this row is offset two spaces, given that each hyperspectral pixel is then stored as `x-value`, `y-value`, and `spectral intensity-values` where the latter has the same length as the list of `Raman Shift` values, created pairs. An exampe .txt fie is provided in `example/example_data.txt` to illustrate this.

The .txt file should be imported into a pandas dataframe.

Although optional, setting the coordinates of the map to start at (0,0) will avoid complications.

```
import pandas as pd

file = 'example/example_data.txt'

file = file.rsplit('.',1)[0]
_df = pd.read_csv(file+'.txt', sep='\t', index_col=False, decimal=',',header=None)
_df = _df.replace(',','', regex=True)
# set origin to zero
_df[0] = _df[0].sub(_df[0].min())
_df[1] = _df[1].sub(_df[1].min())
```

The python modules can then be run in order, for example: 

```
''' --- background reduction --- '''
from raman_hypmap import bgr_process
_df_bgr = bgr_process(_df,file)
_df_bgr = pd.read_csv(file+'_bgr.csv', sep=',', index_col=False,header=None)
_df_raw = pd.read_csv(file+'_raw.csv', sep=',', index_col=False,header=None)

''' --- mineral filter --- '''
from raman_hypmap import filter_window
window_location = '_window_ranges.csv'
filter_window(_df_bgr,file+'_bgr',window_location)

''' --- interactive hyperspectral map --- '''
from raman_hypmap import hypmap_plot
_df_filtered = pd.read_csv(file+'_bgr_w_filter.csv', sep=',', index_col=False)
hypmap_plot(_df_filtered,_df_raw,file)
```

Note that each module imports the generated .csv file from the previous module as an input. These files will be created automatically when running the module to the local folder. 

An example script showing this code is also available in: `example/example.py`

### bgr_process.py

The input format is the standard .txt file generated by the Horiba [LabSpec 6](https://www.horiba.com/int/scientific/products/detail/action/show/Product/labspec-6-spectroscopy-suite-software-1843/) software. However, this format can be recreated via a tab-delimited .txt file as shown below:

| :---         | :---          | :---                                        |
| *blank*      | *blank*       | Raman-shift values (tab-delimited list)     |
| x-coorindate | y-coordinate  | Raman-intensity values (tab-delimited list) |
| :---:        | :---:         | :---:                                       |
| ...          | ...           | ...                                         |

Each spectra is background reduced via an iterative polynomial smoothing method Wang and Dai (2016) and written into a pandas dataframe. Three files are output and automatically saved as .csv files: 
* example_data_raw.csv — reformatted data input saved as a .csv file
* example_data_bg.csv — the background subtracted from each spectra saved as a spectral file
* example_data_bgr.csv — the resulting background reduced spectra saved as a spectral file

Each file is saved in a identical structure to the table shown above, but comma-delimited in a format that can be opened and viewed in Microsoft Excel. 

### filter_window.py

This module needs a .csv file formatted as explained in the section above to operate correctly. The code will itterate over every spectra provided and take the maximum intensity value from this window. It will then output a hyperspectral grid of values using this window. That hyperspecral grid will be saved as an image .png file with greyscale values stretched to the maximum and minimum intensity values recorded. 

The filter windows are imported from a standalone .csv file called `_window_ranges.csv` that can be open and edited in Microsoft Excel. The headings must not be changed. The format of the spredsheet is shown below: 

| phase    | char_min | char_max | comp_min | comp_max | i_min  | i_max  | 
| -------- | -------- |  ------- |  ------- |  ------- |  ----- |  ----- |
| graphite | 1510     | 1620     | 0        | 0        | 0      | 2000   |
| apatite  | 950      | 975      | 0        | 0        | 0      | 2000   |
| :---:    | :---:    | :---:    | :---:    | :---:    | :---:  | :---:  |
| ...      | ...      | ...      | ...      | ...      | ...    | ...    |

These numbers are important as they define how the python code operates:
>'Phase' defines the name of the species
>'char_min' and 'char_max' define the window that the filter will operate over. A new hyperspectral map will be created for every row of values. 
>Sometimes spectra overlap and so 'comp_min' and 'comp_max' define a comparison window. the maximum intensity values defined in the 'char_' window must be *larger* than the maximum intensity value in the 'comp_' window otherwise a value of 0 will be given for that spectra. 
>'i_min' and 'i_max' work in a similar way to the above. The maximum intensity of the spectra in the window defined by 'char_' must be within the intensity thresholds defined by the 'i_min' and 'i_max' variables. 

Each row of values will create and save a new hyperspectral map using that window. 

### hypmap_plot.py



