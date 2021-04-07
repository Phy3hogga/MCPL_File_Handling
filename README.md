# MCPL_File_Handling

Matlab scripts for using MCPL files. Features include
* Reading raw MCPL data into a MAT file format.
* Filtering MAT-formatted MCPL data.
* Writing MCPL data stored in a MAT file into a raw MCPL file format.

### Installation
#### If using HTTPS and not SSH to clone
This repository has submodules linked using SSH rather than HTTPS. Three possible options for cloning the submodules if cloning using HTTPS rather than SSH are as follows:
1. If not already configured, edit the git config file on your PC to re-write all SSH pulls to HTTPS at runtime. Cloning any module with an SSH link should then be automatically redirected to HTTPS.
```git
git config --global url."https://github.com/".insteadOf git@github.com:
git config --global url."https://".insteadOf git://
```
2. Edit the gitmodules file located in the parent repository to manually replace any SSH links with HTTPS link formats. Be aware that when pulling future updates, the gitmodules file may need changing again to update the respective submodules. 
3. Clone each of the the linked submodules manually by opening the individual repositories on Github and cloning them each to their respective sub-directories.

#### Cloning the repository
1. Clone the parent repository MCPL_File_Handling as normal, the submodules will then need to be initiated and cloned seperately. *If using HTTPS to clone this repository, see the above section [If using HTTPS and not SSH](https://github.com/imaging-science-group-ntu/MCPL_File_Handling/readme.md "If using HTTPS and not SSH to clone") before continuing initiating and cloning the required submmodules.*
2. Using the GIT command line from within the parent repository directory, initate all of the submodules using the command.
```git
git submodule update --init --recursive
```
3. Clone all of the submodule contents.
```git
git pull --recurse-submodules
```

#### Compressed MCPL Files
.MCPL files are automatically compressed into a G-Zip format on creation, for automatic unpacking of the <filename>.MCPL.GZ file format it is strongly advised to have WinRAR 5.0 or later installed. In the event that the WinRAR executable (*WinRAR.exe*) is not located on the system enviroment path and fails to be automatically identified as "WinRar.exe", edit the WinRAR_Path variable to point to the appropriate executable. For more information on configuring the WinRAR integration, see [Winrar](https://github.com/Phy3hogga/WinRAR) for a list of valid  optional arguments.
```matlab
%% Parameters for WinRAR implementation
% Path to WinRAR executable (If not automatically found as WinRAR.exe)
RAR_Parameters.WinRAR_Path = 'C:\Program Files\WinRAR\WinRAR.exe';
% Overwrite any files already existing automatically (will not prompt)
RAR_Parameters.Overwrite_Mode = true;
```

### Functions
#### MCPL_To_Mat.m
Converts a binary .MCPL file to a matlab-friendly .MAT file format.
```matlab

```
#### Filter_MCPL_MAT_Data.m
Filters a MAT formatted MCPL file to remove events that are outside a undesired parameter range. Can be useful to reduce insignificant data prior to feeding into another set of simulations.
```matlab

```

#### MAT_To_MCPL.m
Converts a .MAT file format into a binary .MCPL file format.
```matlab

```

#### EKinDir_Pack.m / EKinDir_Unpack.m
performs the compression / decompression algorithms condensing the energy and direction vectors (Dx, Dy, Dz) into three vectors EKinDir_1, EKinDir_2, EKinDir_3 to reduce size on disk. These functions shouldn't need to be directly called, they are used within the scripts listed above.

## Structure of MAT file containing MCPL data
Event data containing multiple arrays (single / double) datatypes where a single integer index between 1 and Header.Particles corresponds to a singular event. Each event has data located at an identical index position from each variable listed below. 
* Position (X, Y, Z)
* Direction Vector (Dx, Dy, Dz)
* Polarisation (Px, Py, Pz) [Optional]
* Energy
* Time
* EKinDir [Packed energy and direction vectors from the MCPL format specification]
* Weight [Optional]
* PDGCode [Optional]
* Userflag [Optional]
* Header [Structure]

Header contains information regarding the MCPL header / data format including
* Endian [Struct] - File endian-ness
* Endian_Switch [Boolean]- If the endian-ness of data needs switching when compared to the native system endian-ness
* MCPL_Version [Integer] - File format of MCPL data, write support is for MCPL version 3, read support version 2 and 3.
* Particles - Total number of events (provides length of each data field)
* Opt_UserFlag [Boolean] - If the optional userflag is being used.
* Opt_Polarisation [Boolean] - If optional polarisation data exists.
* Opt_SinglePrecision [Boolean] - If the data format is single / double.
* Opt_UniversalPDGCode [Int] - If the data only contains a single particle type.
* Opt_ParticleSize [Int] - Number of bytes that contains a single event in data.
* Opt_UniversalWeight [Boolean] - If all events have an identical weighting.
* Opt_Signature [Int] - Currently un-used, part of the C implementation of the MCPL opening processess.
* Source [String] - Source of the data within the MCPL file.
* Comments [String] - Comments
* Blobs [String] - Blobs
* Valid [Boolean] - If the header was read successfully and found to be valid (stops re-checking later).
* End [Int] - Byte position within the MCPL file where the header ends.
* Byte_Size [Int] - Number of bytes that contain one of the dynamic sized fields for an event.
* Byte_Type [String] - Data type used within MCPL data 'single' or 'double'.
* Photon_Byte_Count [Int] - Number of bytes that contain all data related to a single event.
* Byte_Split [Structure] - Contains the relative byte locations for each individual variable relating to a single event in data.
* Sort_Events_By_Weight [Boolean] - On opening an MCPL file, if sorting all events into descending order by most signficant weight in the MAT File.
* Remove_Zero_Weights [Boolean] - If removing all events from the MAT file where event weightings are equal to 0.0.
* File Chunks [Structure] - Used for processing the MCPL file into the MAT file; done in chunks in parallel for multi-core / memory optimisations. 

## Built With

* [Matlab R2018A](https://www.mathworks.com/products/matlab.html)
* [Windows 10](https://www.microsoft.com/en-gb/software-download/windows10)

## References
* ** T. Kittelmann, et al.** - Monte Carlo Particle Lists (MCPL) - [Computer Physics Communications, Volume 218, September 2017, Pages 17-42, ISSN 0010-4655](https://doi.org/10.1016/j.cpc.2017.04.012)
* **Alex Hogg** - *Matlab integration of MCPL* - [Phy3hogga](https://github.com/Phy3hogga)