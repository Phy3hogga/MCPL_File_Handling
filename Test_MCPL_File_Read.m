clear all;
close all;

%Include directories to external GIT libraries
Include_Subdirectories({'Data_Operations','File_Operations','Input_Validation','Parpool', 'Waitbar', 'WinRAR', 'MCPL_Functions'});
%% Test Data
%32 bit data for testing
File_Path = 'F:\MCPL_Output_Diffraction_Test_20210304_145727\MCPL_Output_No_Polarisation_Single_Precision.mcpl.gz';
%64 bit data for testing
%File_Path = 'F:\MCPL_Output_Diffraction_Test_20210303_211044\MCPL_Output_Diffraction_Test_1.mcpl.gz';
%64 bit dataset (small file)
File_Path = 'F:\MCPL_Output_Diffraction_Test_20210329_171051\MCPL_Output_Diffraction_Test_DBL.mcpl.gz';
%File_Path = 'F:\MCPL_Data\XBD_Data\sample_C7H5N3O6_4B_7.xbd';
%File_Path = 'F:\MCPL_Data\Test';

%Prop_Z0
%File_Path = 'F:\MCPL_Output_Diffraction_Test_20210401_143954\MCPL_Monitor_Diffraction_Test_SGL.mcpl.gz';
%No PROP_Z0
%File_Path = 'F:\MCPL_Output_Diffraction_Test_20210401_143359\MCPL_Monitor_Diffraction_Test_SGL.mcpl.gz';


%% Parameters for WinRAR implementation
%Path to WinRAR executable
%RAR_Parameters.WinRAR_Path = 'C:\Program Files\WinRAR\WinRAR.exe';
%By default overwrite any files already existing
RAR_Parameters.Overwrite_Mode = true;

%% Parameters for MCPL processing to MAT file
%If events are sorted in descending order of weight with the most significant events at the top of the file. (true = sort)
Read_Parameters.Sort_Events_By_Weight = true;
%If events with exactly 0 w5eighting (represent no photons) are to be removed (true = removed)
Read_Parameters.Remove_Zero_Weights = true;
%If retaining EKinDir data (reccomended setting to true if wanting to later use subsequent data in simulations)
Read_Parameters.Save_EKinDir = true;
%If the temporary files created during processing are deleted (true = delete temp files)
Read_Parameters.Remove_Temp_Files = false;
%If the GZ archive has already been uncompressed or not (if problems with WinRAR, can bypass decompression)
Read_Parameters.Skip_Uncompress = true;
%Number of cores for the Parpool to use when converting the raw MCPL file (integer)
Read_Parameters.Parpool_Num_Cores = 4;
%Temporary directory to use for constructing / operating on datastore
Read_Parameters.Temp_Directory = 'F:\Windows_Temp_Files';
%Add RAR Parameters to the Read Parameters
Read_Parameters.RAR_Parameters = RAR_Parameters;

Display_Write_Progress = true;
%% Convert MCPL file to MAT file format
tic
Mat_File_Path = MCPL_To_MAT(File_Path, Read_Parameters);
toc
%Perform work on the MAT file
for Current_Mat_File = 1:length(Mat_File_Path)
    %Create new filepath for the MAT file to end up in
    [Directory, Filename, Extension] = fileparts(Mat_File_Path{Current_Mat_File});
    Parent_Directory = fileparts(Directory);
    Output_Directory = fullfile(Parent_Directory, 'DEBUG');
    Attempt_Directory_Creation(Output_Directory);
    MCPL_Filepath = fullfile(Output_Directory, strcat(Filename, '.MCPL'));
    
    %% Filter data within the MAT file
    %Position the events land on the detection plane
    Filters.X.Min = -0.04;
    Filters.X.Max = 0.04;
    Filters.Y.Min = -0.04;
    Filters.Y.Max = 0.04;
    %Filters.Z.Min = 0;
    %Filters.Z.Max = 1;
    %Angle from the normal (to Z), must be +ve valued
    Filters.Angle.Min = 0;
    Filters.Angle.Max = 45;
    %Energy [KeV]
    Filters.Energy.Min = 0;
    Filters.Energy.Max = 130;
    %Weighting
    Filters.Weight.Min = 0.1;
    %Filters.Weight.Max = 35;
    %Create filepath for the filtered MAT file
    Filtered_Mat_File_Path{Current_Mat_File} = fullfile(Output_Directory, strcat(Filename, '-Filtered', Extension));
    
    %Filter the file according to parameters previously set
    tic
    Filtered_Mat_File = MCPL_Filter_MAT_Data(Mat_File_Path{Current_Mat_File}, Filtered_Mat_File_Path{Current_Mat_File}, Filters);
    toc
    
    %Convert MAT file back to an MCPL file
    tic
    MCPL_File = MAT_To_MCPL(Filtered_Mat_File, MCPL_Filepath, Display_Write_Progress);
    toc
    %Test extraction of MCPL file (only for comparison if not filtering data)
    %Mat_File_Path_2 = MCPL_To_MAT(MCPL_File, Read_Parameters);
    % Compare initial MAT and recreated MAT files match
    %visdiff(Mat_File_Path{Current_Mat_File}, Mat_File_Path_2{1});
end