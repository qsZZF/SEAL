Newtopoplot plugin
    A topograhic map visualization plugin for EEGLAB. It is modified from STYLIZED TOPOGRAPHIC MAP PLUGIN
(https://education.msu.edu/kin/hbcl/software.html). The plugin is better compatible with EEGLAB channel location file
compared with the original version and is able to visualize connectivity between channel pairs. A EEGLAB GUI is also 
developed for user-defined plotting.    
============================================================================
Usage:
    GUI Inputs (only part of them are illustrated below): 

        Topo Setting: basic settings for topograph (cartoon head and channels) 
	(check box) Numbers:    plot channel location and their numbers
	(check box) Labels: 	plot channel location and their labels（Only one of Numbers/Labels will work if both selected）
	(check box) Show Inside Only:    plot channels inside the cartoon head only
	(check box) No Electrode:    plot cartoon head only (no electrode would be show)
	(edit) Channel To Plot: default is all

        Topo Weights: settings for spatial weights/power spectrum distribution...
	(edit) EEG Time Range:    plot topograph for EEG signals at specific time points (e.g., EEGdata(: , 20:50))
		- the input could be time point (e.g., 10) or a list of time points (e.g., 10:20, 10 12 13 15)
	(edit) ICA Weights:    plot topographic properties of selected IC(s)
		- e.g., to plot IC 1-5, the input could be 1:5/1 2 3 4 5/1,2,3,4,5
	(edit) Other Data:    user-defined data. The GUI reads the data from basic work space.
		- To use this, first define your data at the work space, e.g., topodata, then the input should be 
		- the variable name "topodata" (without " "). The user-defined data should be a vector/matrix 
		- of [Nchan*1]/[Nchan*m], where Nchan should match the number of channel locations to plot  
	(edit) Color Limit:    default is [0,1]. Note that by default (if you don't edit it) the Data-To-Plot (e.g., EEG Time Range)
		- will be normalized in to [0,1] for each column.
	(check box) Subplot:    works when there are more than one Topo Weights to plot (for data of Nchan*m, m>1). 
		- If selected, they will be shown in one figure using function subplot(). Otherwise, m figures will show.

        Connectivity: settings for connectivity map.
	(edit) Line Width:    used for connectivity lines between channel pairs 
	(edit) Line Color:     likewise,...
	(edit) Other Data:    user-defined data. The GUI reads the data from basic work space (see above).
		- the data should be a matrix/tensor of [Nc*2*Ns]/[Nc*3*Ns], where Nc is the number of connectivity lines.
		- The first two columns define the channel pairs, e.g., [2,5] means a connectivity line (CL) will be drawn 
		- between channel 2 and 5. The third column define the strength of connectivity. Ns ( = 1,2,3...)  defines 
		- the number of connectivity map. 
	(check box) Subplot:     see above
	(check box) Strength:    show strength of connectivity. If selected, the color of CL will be used to reflect the value of
		- strength. In such case, the data must be of [Nc*3*Ns] the (edit) Line Color is not used.  Then the 
		- (edit) Color Limit & (edit) Color Bar could be used as described above.
	(check box) Direction:    show the direction of CL, starts from the first column channel(s) and points to the second 
		- column channel(s).

        others:
	(check box) Plot Topo Weights and Connect In One Fig:
		- If selected, they will be plotted in on fig. If more than one topo weights & connectivity maps are to plot,
		- the (check box) Subplot / (checkbox) Show Inside Only in Topo Weights settings/ Connectivity settings
		- work the same way as described above.  Note that in this case the (check box) Subplot is shared by 
		- Topo Weights and Connect, which means that the Subplot can work if one of the check box is selected.
		- Meanwhile, other settings (e.g., (edit) Color Limit for Topo & Connect) could be defined with different 
		- parameters.
	(check box) Plot Gif:
		- not used yet.

	**See "help Newtopoplot" for detailed information**
	The users are also recommend to call the function Newtopoplot from command line
=================================================================================
    The Newtopoplot plugin has not been fully developed yet and should only be used within the lab for now.      
