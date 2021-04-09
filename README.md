# MCPL_File_Handling

A series of Matlab scripts for translating and manipulating **M**onte-**C**arlo **P**article **L**ists (**MCPL**) files within Matlab.

<table>
	<thead>
		<tr>
			<th>Script</th>
			<th>Purpose</th>
			<th>Example</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td>MCPL_To_MAT.m</td>
			<td>Translating binary .MCPL data into a .MAT file format.</td>
			<td>
				<a href="README.md#mcpl_to_mat" title="Filter">Example</a>
			</td>
		</tr>
		<tr>
			<td>Filter_MCPL_MAT_Data.m</td>
			<td>Creates a filtered copy of a MAT file containing MCPL data within specific constraints</td>
			<td>
				<a href="README.md#Filter_MCPL_MAT_Data" title="Filter">Example</a>
			</td>
		</tr>
		<tr>
			<td>MAT_To_MCPL.m</td>
			<td>Translating a .MAT file format into a binary .MCPL file format</td>
			<td>
				<a href="README.md#mat_to_mcpl" title="Filter">Example</a>
			</td>
		</tr>
	</tbody>
</table>

## Functions
### MCPL_To_Mat
Converts a binary .MCPL file to a matlab-friendly .MAT file format.

***Note:** WinRAR implementation is currently only supported on windows. On other operating systems compressed .MCPL.GZ files will need to be uncompressed manually and Skip_Uncompress will need setting to true.*
```matlab
%% File path to convert the MCPL file to MAT format
File_Path = 'D:\MCPL_Monitor_Diffraction_Test_SGL.mcpl.gz';

%% Parameters for WinRAR implementation
% Path to WinRAR executable
RAR_Parameters.WinRAR_Path = 'C:\Program Files\WinRAR\WinRAR.exe';
% By default overwrite any files already existing
RAR_Parameters.Overwrite_Mode = true;

%% Parameters for MCPL processing to MAT file
% If events are sorted in order of weighting with the most significant events at the top of the file. (true = sort)
Read_Parameters.Sort_Events_By_Weight = true;
% If events with exactly 0 weighting (represent no photons) are to be removed (true = removed)
Read_Parameters.Remove_Zero_Weights = true;
% Number of cores for the Parpool to use when converting the raw MCPL file (integer)
Read_Parameters.Parpool_Num_Cores = 6;
% If the temporary files created during multi-core processing are deleted (true = delete temp files)
Read_Parameters.Remove_Temp_Files = true;
% Add RAR Parameters to the Read Parameters to pass it through to decompress the GZip archive
Read_Parameters.RAR_Parameters = RAR_Parameters;
% If the GZ archive has already been uncompressed, skip decompression.
% If problems with WinRAR occur, can bypass decompression (true = disable decompression)
Read_Parameters.Skip_Uncompress = false;

%% Convert MCPL file to MAT file format.
% Mat_File_Path returns a cell array containing each individual file within the GZIP archive.
% This compensates for the case of non-merged MPI data).
Mat_File_Path = MCPL_To_MAT(File_Path, Read_Parameters);
```
### Filter_MCPL_MAT_Data
Filters a MAT formatted MCPL file to remove events that are outside a undesired parameter range. Can be useful to reduce insignificant data prior to feeding into another set of simulations. Any filters that aren't assigned will not be applied.
```matlab
%% For each MAT file created by MCPL_To_MAT
for Current_Mat_File = 1:length(Mat_File_Path)
	%% Filtering parameters to discard data within the MAT file
	% Position that events land on the detection plane
	Filters.X.Min = -0.04;
	Filters.X.Max = 0.04;
	Filters.Y.Min = -0.04;
	Filters.Y.Max = 0.04;
	Filters.Z.Min = 0;
	Filters.Z.Max = 1;
	% Angle from the normal (to Z) [positive valued]
	Filters.Angle.Min = 10;
	Filters.Angle.Max = 45;
	% Energy [KeV]
	Filters.Energy.Min = 5;
	Filters.Energy.Max = 130;
	% Weighting
	Filters.Weight.Min = 0.05;
	Filters.Weight.Max = 35;

	%% Create filepath to save the filtered MAT file
	Filtered_Mat_File_Path{Current_Mat_File} = strcat(Mat_File_Path{Current_Mat_File}, '-Filtered');

	%% Filter the MCPL data within the MAT file
	Filtered_Mat_File = Filter_MPCL_MAT_Data(Mat_File_Path{Current_Mat_File}, Filtered_Mat_File_Path{Current_Mat_File}, Filters);
end
```

### MAT_To_MCPL
Converts a .MAT file format into a binary .MCPL file format.
```matlab
%% For each MAT file created by MCPL_To_MAT
for Current_Mat_File = 1:length(Mat_File_Path)
	%% Create directory path to save the .MAT file to
	[Directory, Filename, Extension] = fileparts(Mat_File_Path{Current_Mat_File});
	Parent_Directory = fileparts(Directory);
	% Create new directory path to save converted file to
	Output_Directory = fullfile(Parent_Directory, 'Converted');
	Attempt_Directory_Creation(Output_Directory);
	% Create file path
	MCPL_Filepath = fullfile(Output_Directory, strcat(Filename, '.MCPL'));

	%% Convert MAT file back to an MCPL file
	MCPL_File = MAT_To_MCPL(Mat_File_Path{Current_Mat_File}, MCPL_Filepath);
end
```

### EKinDir_Pack / EKinDir_Unpack
These scripts perform the compression / decompression algorithms condensing the energy and direction vectors (Dx, Dy, Dz) into three vectors EKinDir_1, EKinDir_2, EKinDir_3 to reduce filesizes. These functions shouldn't need to be directly called, as they are exclusively used within the scripts listed above when appropriate.

## Structure of MAT file containing MCPL data
### Event Data
Event data containing multiple individual arrays which are all universally indexed between 1 and Header.Particles, where an individual index corresponds to a singular event. The data present in the MAT file is shown in the table directly below, the header structure is shown in the subsequent table. 

<table>
	<thead>
		<tr>
			<th>Variable</th>
			<th>Label</th>
			<th>Datatype</th>
			<th>Optional</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td>Header</td>
			<td>MCPL File Header Data</td>
			<td>Structure</td>
			<td>No (See <a href="README.md#Header" title="Header Table">header table</a>)</td>
		</tr>
		<tr>
			<td>X</td>
			<td rowspan=3>Position [m]</td>
			<td rowspan=3>Single/Double*</td>
			<td rowspan=3>No</td>
		</tr>
		<tr>
			<td>Y</td>
		</tr>
		<tr>
			<td>Z</td>
		</tr>
		<tr>
			<td>Dx</td>
			<td rowspan=3>Direction Vector</td>
			<td rowspan=3>Single/Double*</td>
			<td rowspan=3>No</td>
		</tr>
		<tr>
			<td>Dy</td>
		</tr>
		<tr>
			<td>Dz</td>
		</tr>
		<tr>
			<td>Px</td>
			<td rowspan=3>Polarisation</td>
			<td rowspan=3>Single/Double*</td>
			<td rowspan=3>Yes</td>
		</tr>
		<tr>
			<td>Py</td>
		</tr>
		<tr>
			<td>Pz</td>
		</tr>
		<tr>
			<td>Energy</td>
			<td>Energy [KeV]</td>
			<td>Single/Double*</td>
			<td>No</td>
		</tr>
		<tr>
			<td>Time</td>
			<td></td>
			<td>Single/Double*</td>
			<td>No</td>
		</tr>
		<tr>
			<td>Weight</td>
			<td>Event Weighting</td>
			<td>Single/Double*</td>
			<td>Yes (If not universally weighted)</td>
		</tr>
		<tr>
			<td>PDGCode</td>
			<td>Particle idenfitier</td>
			<td>Integer</td>
			<td>Yes (If varying particle types)</td>
		</tr>
		<tr>
			<td>Userflag</td>
			<td>User Specified Flag</td>
			<td>Integer</td>
			<td>Yes (If value set)</td>
		</tr>
		<tr>
			<td>EKinDir_1</td>
			<td rowspan=3>Compressed Energy + Direction Vectors</td>
			<td rowspan=3>Single/Double*</td>
			<td rowspan=3>Yes (Not included by default, can be recalculated from Energy and Direction vectors)</td>
		</tr>
		<tr>
			<td>EKinDir_2</td>
		</tr>
		<tr>
			<td>EKinDir_3</td>
		</tr>
	</tbody>
</table>

\* Data type varies between single or double datatypes due to the precision specified when originally creating the MCPL file.

### Header
Header contains the information regarding either the MCPL header / data format or miscilaneous values used by multiple subsequent functions in this software package. If the Header variable is missing from the MAT file created by MCPL_To_MAT, the MCPL to MAT file conversion was unsuccessful; as the header is the final part of the MAT file to be written.

<table>
	<thead>
		<tr>
			<th>Structure Field</th>
			<th>Datatype</th>
			<th>Notes</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td>Endian</td>
			<td>Struct</td>
			<td>File Endian-ness</td>
		</tr>
		<tr>
			<td>Endian_Switch</td>
			<td>Boolean</td>
			<td>System Endian-ness matches File Endian-ness at the time of translation</td>
		</tr>
		<tr>
			<td>MCPL_Version</td>
			<td>Integer</td>
			<td>Supports uncompressing MCPL versions 2 and 3, compresses version 3</td>
		</tr>
		<tr>
			<td>Particles</td>
			<td>Integer</td>
			<td>Total number of events</td>
		</tr>
		<tr>
			<td>Opt_Userflag</td>
			<td>Boolean</td>
			<td>If the optional userflag is set</td>
		</tr>
		<tr>
			<td>Opt_Polarisation</td>
			<td>Boolean</td>
			<td>If polarisation data exists</td> 
		</tr>
		<tr>
			<td>Opt_SinglePrecision</td>
			<td>Boolean</td>
			<td>If single / double precision</td>
		</tr>
		<tr>
			<td>Opt_UniversalPDGCode</td>
			<td>Integer</td>
			<td>PDGCode identifying a single particle type if no variance</td>
		</tr>
		<tr>
			<td>Opt_ParticleSize</td>
			<td>Integer</td>
			<td>Number of bytes that contain a single event</td>
		</tr>
		<tr>
			<td>Opt_UniversalWeight</td>
			<td>Boolean</td>
			<td>If all events have identical weighting</td>
		</tr>
		<tr>
			<td>Opt_Signature</td>
			<td>Int</td>
			<td>Currently un-used, part of the C implementation of the MCPL opening process</td>
		</tr>
		<tr>
			<td>Source</td>
			<td>Cell Array</td>
			<td>Source of the MCPL data</td>
		</tr>
		<tr>
			<td>Comments</td>
			<td>Cell Array</td>
			<td>Comments in the MCPL data</td>
		</tr>
		<tr>
			<td>Blobs</td>
			<td>Cell Array</td>
			<td>Blobs in the MCPL data</td>
		</tr>
		<tr>
			<td>Valid</td>
			<td>Boolean</td>
			<td>If the header was read successfully and found to be a valid format</td>
		</tr>
		<tr>
			<td>End</td>
			<td>Integer</td>
			<td>Final Byte position of the header data</td>
		</tr>
		<tr>
			<td>Byte_Size</td>
			<td>Integer</td>
			<td>Number of bytes that contain a single dynamicly sized variable</td>
		</tr>
		<tr>
			<td>Byte_Type</td>
			<td>String</td>
			<td>Data type used within the dynamic size MCPL files ('Single' / 'Double')</td>
		</tr>
		<tr>
			<td>Photon_Byte_Count</td>
			<td>Integer</td>
			<td>Number of bytes that contain all data related to a single event</td>
		</tr>
		<tr>
			<td>Byte_Split</td>
			<td>Structure</td>
			<td>Contains the relative byte positions for all variables corresponding to a single event</td>
		</tr>
		<tr>
			<td>Sort_Events_By_Weight</td>
			<td>Boolean</td>
			<td>On opening an MCPL file, if sorting all events into descending order by most signficant weight in the MAT File.</td>
		</tr>
		<tr>
			<td>Remove_Zero_Weights</td>
			<td>Boolean</td>
			<td>If removing all insignificant events from the MAT file where event weightings are equal to exactly 0.</td>
		</tr>
		<tr>
			<td>File_Chunks</td>
			<td>Structure</td>
			<td>Tracking tempoary files, byte positions in the original file, number of events</td>
		</tr>
	</tbody>
</table>


## Installation Requirements / Advice
### Compressed MCPL Files (.MCPL.GZ)
***Note:** WinRAR implementation is currently only supported on windows.*

.MCPL files are automatically compressed into a G-Zip format on creation, for automatic unpacking of the <filename>.MCPL.GZ file format it is strongly advised to have WinRAR 5.0 (or later) installed. In the event that the WinRAR executable (*WinRAR.exe*) is not located on the system enviroment path and fails to be automatically identified as "WinRar.exe", edit the WinRAR_Path variable to point to the appropriate executable. For more information on configuring the WinRAR integration, see the [WinRAR submodule readme](https://github.com/Phy3hogga/WinRAR) for a list of addditional optional arguments.
```matlab
%% Parameters for WinRAR implementation
% Path to WinRAR executable (If not automatically found as WinRAR.exe)
RAR_Parameters.WinRAR_Path = 'C:\Program Files\WinRAR\WinRAR.exe';
% Overwrite any files already existing automatically (will not prompt)
RAR_Parameters.Overwrite_Mode = true;
%Add RAR Parameters to the Read Parameters for MPCL_To_MAT.m
Read_Parameters.RAR_Parameters = RAR_Parameters;
```

### Cloning the repository
1. Clone the parent repository MCPL_File_Handling as normal, the submodules will then need to be initiated and cloned seperately. *If using HTTPS to clone repositories from github, see section [If using HTTPS and not SSH](README.md#if-using-https-and-not-ssh-to-clone) before attempting to initiate and clone the required submmodules.*
2. Using the GIT command line from within the parent repository directory, initate all of the submodules using the command. This command only needs performing before the pulling the submodule contents for the first time.
```git
git submodule update --init --recursive
```
3. Clone all of the submodule contents.
```git
git pull --recurse-submodules
```

### If using HTTPS and not SSH to clone
This repository has submodules linked using SSH rather than HTTPS. If wanting to clone the repository using HTTPS, additional steps to clone the submodules using HTTPS rather than SSH are as follows:
* If not already configured, edit the git config file on your PC to re-write all SSH pulls to HTTPS at runtime. Cloning any module with an SSH link should then be automatically redirected to HTTPS.
```git
git config --global url."https://github.com/".insteadOf git@github.com:
git config --global url."https://".insteadOf git://
```
* Edit the gitmodules file located in the parent repository to manually replace any SSH links with HTTPS link formats. Be aware that when pulling future updates, the gitmodules file may need changing again to update the respective submodules.
* Clone each of the the linked submodules manually by opening the individual repositories on Github and cloning them each to their respective sub-directories.

## Built With

* [Matlab R2018A](https://www.mathworks.com/products/matlab.html)
* [WinRAR 5.0+](https://www.rarlab.com/)
* [Windows 10](https://www.microsoft.com/en-gb/software-download/windows10)

## References
* **Alex Hogg** - *Matlab integration of MCPL* - [Phy3hogga](https://github.com/Phy3hogga)
* **T. Kittelmann, et al.** - *Monte Carlo Particle Lists (MCPL)* - [Computer Physics Communications, Volume 218, September 2017, Pages 17-42, ISSN 0010-4655](https://doi.org/10.1016/j.cpc.2017.04.012)