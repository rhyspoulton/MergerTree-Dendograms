# Mergertree Dendogram Builder

Developed by:
Rhys Poulton

This is a code to create the MergerTree Dendograms as presented in [Poulton et at. 2018](https://arxiv.org/abs/1809.06043).

## Setup

The code uses a submodule so for a fresh clone use:

`git submodule update --init --recursive`

to get the submodule or to update it use:

`git submodule update --recursive --remote`

## Running

First the halo catalogue needs to be converted into Efficient Tree Format (ETF) to be readable by the dendogram code, there exist conversion tools for VELOCIraptor, AHF and Rockstar. If the halo catalogue is in one of these formats then please edit the config file convToETF.cfg supplied in the directory so the necessary values are set. To load in the data from the different halo catalogues the requires different information depending on the catalogue:

- **VELOCIraptor** - For this halo catalogue the code requires the unified halo catalogue produced from the python tools
- **AHF** - Requires a file containing a list of the AHF \*halo output filenames and another file with the list of the MergerTree \*_mtree filenames. These lists need to be ascending in redshift.

- **Rockstar** - A filelist containing the Consistent Trees *.dat filenames

After the config file has been altered, the conversion tool can be run by typing:

`python convToETF.py <Format> convToETF.cfg`

Where \<Format> is the merger tree format which can be either AHF, Rock (Rockstar) or VEL (VELOCIraptor). There is a additional converter build for Millennium (Mill) format, but due to the slightly different formats used by different people this may need to be changed. The name of the file that the converted data will be put in with the name \<outfilename>.ETF.tree.hdf, given in the convToETF.cfg.

Once the file has been converted then the dendograms can be built, this can be done by typing:

`python CreateDendogram.py <ETFFilename> <nplot> <outputdir> Plot_config.cfg`

Where \<ETFFilename> is the name for the converted ETF file, \<nplot> is the number of dendograms to be plotted, \<outputdir> is the output folder to put the dendograms (will be created if does not exit) and the plot_config.cfg is the config file containing all the options for plotting. The code will ask what dataset to set the size of the points and if plotColorBar = 1 is set then it will also ask what dataset to set as the color. If this is zero then the code will set the color based on if the halo is a subhalo or a halo. If overplotdata = 1 the code will also ask what data to set as the data over the plot for each halo.

Once that has run you should have your first \<nplot> dendograms in the \<outputdir> folder.

## Notes

This code is still currently under development so please contact me if you run into any problems or you think something can be improved