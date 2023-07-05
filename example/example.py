import pandas as pd
import os
import locale
locale.setlocale(locale.LC_NUMERIC, "en_DK.UTF-8")

''' --- import file --- '''
file = 'example/example_data.txt'
file = file.rsplit('.',1)[0]
# import RAW .txt file into dataframe
_df = pd.read_csv(file+'.txt', sep='\t', index_col=False, decimal=',',header=None)
_df = _df.replace(',','', regex=True)
# set origin to zero
_df[0] = _df[0].sub(_df[0].min())
_df[1] = _df[1].sub(_df[1].min())

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

exit()


