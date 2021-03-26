clear all;
close all;
Include_Subdirectories({'Parpool','WinRAR','File_Operations','Geometric'});
File_Path = 'C:\Users\alex_\Desktop\MCPL_Output_Diffraction_Test_20210302_123431\MCPL_Output_Diffraction_Test.mcpl.gz';
%32 bit
File_Path = 'D:\MCPL_Output_Diffraction_Test_20210304_145727\MCPL_Output_No_Polarisation_Single_Precision.mcpl.gz';
%64 bit
%File_Path = 'D:\MCPL_Output_Diffraction_Test_20210303_211044\MCPL_Output_Diffraction_Test_1.mcpl.gz';

%% Read parameters

Read_Parameters.Sort_Events_By_Weight = true;
Read_Parameters.Remove_Zero_Weights = true;
Read_Parameters.Remove_Temp_Files = true;
Read_Parameters.Skip_Uncompress = false;
Read_Parameters.Parpool_Num_Cores = 6;

%Convert binary MCPL file to MAT file format
Mat_File_Path = MCPL_To_MAT(File_Path, Read_Parameters);