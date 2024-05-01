import scipy.io
import pandas as pd

# This script takes a .mat file created by the calcium_imaging.m script and turns it into a binary pandas dataframe so it can be analyzed by python scripts
# You may want to iterate through multiple output .mat files if you started with more than one .tif file


mat = scipy.io.loadmat('test001.mat') #rename to whatever your mat file is called
spikeTrain = mat['data'][0][0][0]


spikeTrain_binary = (spikeTrain != 0).astype(int)
spikeTrain_df = pd.DataFrame(spikeTrain_binary)


print(spikeTrain_df)
