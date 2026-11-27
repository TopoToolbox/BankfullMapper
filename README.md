[![BankfullMapper-for-MATLAB on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/181632-bankfullmapper)

# BankfullMapper

<img src="images/bankfull_map.png" align="center" width=500>

BankfullMapper is a semi-automatic toolkit for bankfull geometry extraction and discharge estimation.


## Installation

You can either download the entire toolbox and add the folder `toolbox` and its subfolders to the search path. Alternatively, download the latest release of the toolbox (.mltbx). Once, you double-click the toolbox, it will be installed to MATLAB's add-on folder.

This toolkit includes MATLAB functions to analyze the river bankfull geometry and derive discharge estimate applying Manning's equation. 

This toolkit requires TopoToolbox, which can be downloaded from [here](https://github.com/TopoToolbox/topotoolbox3), and needs the same MATLAB's built-in packages (Image Processing Toolbox, Mapping Toolbox) required for running TopoToolbox.


Peak detection uses MATLAB's `islocalmax`, so the Signal Processing Toolbox is **not** required; `ksdensity` and `knnsearch` require the Statistics and Machine Learning Toolbox.


## What's new (November 2026)

- **Closed profiles and spill elevation (`bankfull`, `detect_peak`, `manningseq`)**: each profile is closed at both ends with vertical walls higher than `max_depth`, so the water is always confined within the profile. This fixes a bug in which the flow width could be computed on the *dry* portions of the profile (when a profile end was below the water level), producing spurious peaks of the hydraulic depth curve, blocks of flooded floodplain and gaps in the maps. `bankfull` returns the spill elevation of each profile (`bank.zspill`, the lowest of the two highest points on the left and right of the centerline) and warns when `zspill < max_depth`: in that case the profile does not contain the valley sides and its `width` in `prof` should be increased (3-5 times the channel width is a good starting point).
- **`'connected'` option** (default `true`) in `bankfull`, `detect_peak`, `manningseq` and `visual`: only the wet portion of the profile containing the centerline, and in the maps only the flooded cells connected to the river course, are considered. Ditches, ponds, oxbows or the other limb of a meander are excluded. Use `'connected',false` (and `'clip',false` in `prof`, `'nearest',false` in `visual`) to reproduce the original behaviour.
- **Tight meanders**: `prof` has a new `'clip'` option (default `true`) that cuts each transverse profile where it reaches another reach of the river (e.g. the other limb of a tight meander); `visual` has a new `'nearest'` option (default `true`) that assigns each cell to the nearest point of the river course instead of interpolating the along-course distance from the swath, so overlapping transverse profiles no longer produce misassigned cells.
- **`'thalwegwindow'` option** (default 10 m) in `bankfull`, `detect_peak`, `manningseq` and `visual`: the thalweg of each profile is the lowest point within this distance from the river course, and the connected wet area starts from it. Increase it when the river course does not follow the channel, e.g. when the profiles of one river course are used to compare DEMs of different dates after channel migration (see `example_code2.m`). `bankfull` returns the position and elevation of the thalweg (`bank.dthalweg`, `bank.zthalweg`) and warns when the thalweg is not on the river.
- **`visual`**: each profile is applied to a reach centred on the profile (boundaries halfway between consecutive profiles).
- **NaN-safe profiles**: `bankfull`, `detect_peak` and `manningseq` handle profiles with NaN values (clipped ends or nodata cells).
- **Exact and faster cross-section geometry**: flow area, width and wetted perimeter are computed exactly on the vertices of each profile and on its intersections with the water level, instead of interpolating each profile on 10^6 points. The results are exact and the processing is much faster (in particular in `manningseq`).
- **`manningseq`**: fixed the hydraulic radius (element-wise `./` instead of matrix division `/`), which affected all discharge values.
- **`detect_peak`**: no longer uses the Signal Processing Toolbox; fixed errors when two peaks have the same prominence (the lowest is taken) or two levels have the same hydraulic depth.
- **`sectiongradient`**: the flow direction of the river course is detected automatically; the first section no longer has a NaN slope.
- **`har`**: the REM is computed by Inverse Distance Weighting (IDW) of the k nearest river points (optional `'power'`, `'k'`, `'plot'`); the elevation of the river course is sampled from the DEM, so a STREAMobj derived from a different grid can be used.
- **`prof`**: new `'keeptrace'` option (default `true`): when the trace is smoothed, the transverse profiles are oriented by the smoothed trace but remain centred on the original river course (with strong smoothing the profile centre could fall far from the channel). 
- **`prof`**: `'smooth'` is now a switch (`true`/`false`) plus `'smoothingfactor'` (0-1), consistent with `SWATHobj` in TopoToolbox 3.
- **`WORKFLOW.m`**: the slope computed on a STREAMobj is mapped to the sections.

Note: results can differ from those obtained with previous versions (in particular with `'connected',true` and with the new REM).


## Getting Started

In the examples folder, you can find example_code.m and example_code2.m where two examples are reported.

The data for the examples can be downloaded from this [Google Drive](https://drive.google.com/drive/folders/1dM2VZWRmLslaq0dfifNAPaci4E1Bd0Pa?usp=drive_link).

To run the tool with your data, just follow the general workflow in the WORKFLOW.m script. 

## References

Delchiaro, M., Ruscitto, V., Schwanghart, W., Brignone, E., Piacentini, D., and Troiani, F.. “BankfullMapper: a semi-automated MATLAB tool on high-resolution Digital Terrain Models for spatio-temporal monitoring of bankfull geometry and discharge”. Computers and Geosciences. https://doi.org/10.1016/j.cageo.2025.106001

Schwanghart, W., & Scherler, D., 2014. “TopoToolbox 2–MATLAB-based software for topographic analysis and modeling in Earth surface sciences”. Earth Surface Dynamics, 2(1), 1-7. https://doi.org/10.5194/esurf-2-1-2014
