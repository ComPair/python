Various scripts and config files to obtain simulated background spectra from `.tra` files.

`RemoveBottomTrackerLeptonicInteractions.py`: This should be called once after the background sims are run. Removes background component of trapped leptonic events that come in via the gap below the ACD.

`mimrec_bg_spectra.sh`: Runs mimrec over tra files to produce root files with background spectra. Input/output directory, config files and detector geometry are all hardcoded and should be adjusted. This produces one root file per component (TrappedLeptonic, Photonic etc) and event type (Tracked Compton, Pair etc). The spectra are not normalized by exposure time yet!

`PlotBackground.cxx`: Reads in the root files from the previous step to normalize and plot the spectra and calculate the total background rates (sum over all components). Again, input and output files and simulated time are hardcoded. Output: plots (png, pdf) of the reconstructed background spectrum for each component, root files with said spectra, and root files with just the total background histogram for each event type. Use the `hadd` macro to combine those into one file if desired. 

`PlotTotalBG.C`: Root script to plot total Untracked Compton/Tracked Compton/Pair background spectra (summed over all components), using the output from the previous script.

`mimrec_AMEGOX_PlotBkgSpectrum_all.cfg`: Config file, selecting both Compton and pair events.

`mimrec_AMEGOX_PlotBkgSpectrum_C.cfg`: Config file, Compton only.

`mimrec_AMEGOX_PlotBkgSpectrum_P.cfg`: Config file, Pair only.

`mimrec_AMEGOX_PlotBkgSpectrum_TC.cfg`: Config file, tracked Compton only.

`mimrec_AMEGOX_PlotBkgSpectrum_UC.cfg`: Config file, untracked Compton only.

