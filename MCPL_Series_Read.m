clear all;
close all;

%Include directories to external GIT repositories 
Include_Subdirectories({'Data_Operations','File_Operations','Input_Validation','Parpool', 'Waitbar', 'WinRAR', 'MCPL_Functions'});

%% File type for searching for raw data
%File_Type = '.mcpl';
File_Type = '.xbd';

%% Test Data
%32 bit data for testing
Directory_Path = 'D:\Simulations\81\out';
%Directory_Path = 'D:\Data_Capture\0Deg\data';
List_File_Path = Search_Files(Directory_Path, File_Type);
%Output MAT file path
Merged_Filename = 'merged';

[Root_File_Path, Filename, ~] = fileparts(Directory_Path);
Merged_File_Path = fullfile(Root_File_Path, strcat(Merged_Filename, ".mat"));

%% Parameters for WinRAR implementation
%Path to WinRAR executable
%RAR_Parameters.WinRAR_Path = 'C:\Program Files\WinRAR\WinRAR.exe';
%By default overwrite any files already existing
RAR_Parameters.Overwrite_Mode = true;

%% Parameters for MCPL processing to MAT file
%If events are sorted in descending order of weight with the most significant events at the top of the file. (true = sort)
Read_Parameters.Sort_Events_By_Weight = true;
%If events with exactly 0 weighting (represent no photons) are to be removed (true = removed)
Read_Parameters.Remove_Zero_Weights = true;
%If retaining EKinDir data (reccomended setting to true if wanting to later use subsequent data in simulations)
Read_Parameters.Save_EKinDir = true;
%If the GZ archive has already been uncompressed or not (if problems with WinRAR, can bypass decompression)
Read_Parameters.Skip_Uncompress = true;
%Number of cores for the Parpool to use when converting the raw MCPL file (integer)
Read_Parameters.Parpool_Num_Cores = 6;
%Temporary directory to use for constructing / operating on datastore
Read_Parameters.Temp_Directory = 'D:\Windows_Temp_Files';
%Add RAR Parameters to the Read Parameters
Read_Parameters.RAR_Parameters = RAR_Parameters;

Display_Write_Progress = true;
%% Convert MCPL file to MAT file format
Mat_File_Path = MCPL_To_MAT(Directory_Path, Read_Parameters);

%% If multiple files provided by MCPL_To_MAT, merge.
Mat_File_Path = MCPL_Merge_Files(Mat_File_Path, Merged_File_Path, false);

%% Load data
%read(Mat_File_Path);
%Mat_File_Path = "D:\Data_Capture\Alex_Hex2\sample.mat";
%File_Data_Store = tall(fileDatastore(Mat_File_Path, 'ReadFcn', @(x)struct2table(load(x, 'Dx', 'Dy', 'Dz', 'Energy', 'Weight', 'X', 'Y', 'Z')), 'UniformRead', true));
File_Content = load(Mat_File_Path, 'X', 'Y', 'Dz', 'Weight', 'Energy');
X = File_Content.X;
Y = File_Content.Y;
Dz = File_Content.Dz;
Weight = File_Content.Weight;
Energy = File_Content.Energy;

%% Calculations and plots
%Calculate vector from Z normal; simplification of dot product with unit vectors.
Event_Angle = acosd(Dz);
%Event_Angle = gather(Event_Angle);

% %Display X, Y, Z output data from source
% figure();
% scatter3(X, Y, Z, [], Event_Angle, '.');
% xlabel('X [m]');
% ylabel('Y [m]');
% 
% %Display Dz, Dy, Dz output data from source
% figure();
% scatter3(Dx, Dy, Dz, [], Event_Angle, '.');
% xlabel('X [m]');
% ylabel('Y [m]');
% 
% %Display histogram of Angular data from source
% figure();
% histogram(Event_Angle);
% xlabel(['Incident Angle [', char(176), ']']);
% ylabel('Frequency');
% 
% %Display event angle as a 3d surface
% figure();
% scatter3(X, Y, Event_Angle, [], Event_Angle, '.');
% xlabel('X [m]');
% ylabel('Y [m]');

%Rebinning parameters
Num_Spatial_Bins = 150;
Num_Hist_Bins = 150;
Spatial_Bin_Tol = 0.1e-3;
%Weight = gather(File_Data_Store.Weight);
%X = gather(File_Data_Store.X);
%Y = gather(File_Data_Store.Y);

%Create equally spaced bins in X and Y
X_Bins = linspace(min(X(:)) - Spatial_Bin_Tol, max(X(:)) + Spatial_Bin_Tol, Num_Spatial_Bins + 1);
Y_Bins = linspace(min(Y(:)) - Spatial_Bin_Tol, max(Y(:)) + Spatial_Bin_Tol, Num_Spatial_Bins + 1);
%Create grid from bins
[Grid_X, Grid_Y] = ndgrid(X_Bins, Y_Bins);
Weighted_Binned_Angle_Min = nan(size(Grid_X));
Weighted_Binned_Angle_Max = nan(size(Grid_X));
Weighted_Binned_Angle_Mean = nan(size(Grid_X));
Weighted_Binned_Angle_Weighted_Mean = nan(size(Grid_X));
Weighted_Binned_Angle_Std = nan(size(Grid_X));
Weighted_Binned_Angle_Weighted_Std = nan(size(Grid_X));
Weighted_Binned_Angle_Count = nan(size(Grid_X));

Parpool_Create();
for Current_X = 1:length(X_Bins) - 1
    %Get X Index
    X_Index = (X_Bins(Current_X) < X) & (X <= X_Bins(Current_X + 1));
    %Preallocate line data
    Weighted_Binned_Angle_Line_Min = nan(size(Y_Bins));
    Weighted_Binned_Angle_Line_Max = nan(size(Y_Bins));
    Weighted_Binned_Angle_Line_Mean = nan(size(Y_Bins));
    Weighted_Binned_Angle_Line_Std = nan(size(Y_Bins));
    Weighted_Binned_Angle_Line_Weighted_Mean = nan(size(Y_Bins));
    Weighted_Binned_Angle_Line_Weighted_Std = nan(size(Y_Bins));
    Weighted_Binned_Angle_Line_Count = nan(size(Y_Bins));
    parfor Current_Y = 1:length(Y_Bins) - 1
        % Get Y index
        Y_Index = (Y_Bins(Current_Y) < Y) & (Y <= Y_Bins(Current_Y + 1));
        % Combine X and Y index limits
        Index = X_Index & Y_Index;
        % Load data
        Pixel_Non_Weighted_Angle = Event_Angle(Index);
        Pixel_Weight = Weight(Index);
        if(~isempty(Pixel_Non_Weighted_Angle))
            %% Calculations
            %Weighted_Angle = Dz(Index);
            Weighted_Angle = sum((Pixel_Non_Weighted_Angle .* Pixel_Weight), 'omitnan') ./ sum(Pixel_Weight, 'omitnan');
            Number_Nonzero_Weights = sum(Floating_Point_Equal(Pixel_Weight, 0) == 0);
            Weighted_Std = sqrt(abs(sum(Pixel_Weight .* (Pixel_Non_Weighted_Angle - Weighted_Angle).^2,'omitnan')/(((Number_Nonzero_Weights - 1)/Number_Nonzero_Weights) .* sum(Pixel_Weight, 'omitnan'))));
            %% Hold data for single X line
            Weighted_Binned_Angle_Line_Min(Current_Y) = min(Pixel_Non_Weighted_Angle);
            Weighted_Binned_Angle_Line_Max(Current_Y) = max(Pixel_Non_Weighted_Angle);
            Weighted_Binned_Angle_Line_Mean(Current_Y) = mean(Pixel_Non_Weighted_Angle, 'omitnan');
            Weighted_Binned_Angle_Line_Std(Current_Y) = std(Pixel_Non_Weighted_Angle, 'omitnan');
            Weighted_Binned_Angle_Line_Weighted_Mean(Current_Y) = Weighted_Angle;
            Weighted_Binned_Angle_Line_Weighted_Std(Current_Y) = Weighted_Std;
            Weighted_Binned_Angle_Line_Count(Current_Y) = length(Pixel_Non_Weighted_Angle);
        end
    end
    %% Save data
    Weighted_Binned_Angle_Min(Current_X, :) = Weighted_Binned_Angle_Line_Min;
    Weighted_Binned_Angle_Max(Current_X, :) = Weighted_Binned_Angle_Line_Max;
    Weighted_Binned_Angle_Mean(Current_X, :) = Weighted_Binned_Angle_Line_Mean;
    Weighted_Binned_Angle_Std(Current_X, :) = Weighted_Binned_Angle_Line_Std;
    Weighted_Binned_Angle_Weighted_Mean(Current_X, :) = Weighted_Binned_Angle_Line_Weighted_Mean;
    Weighted_Binned_Angle_Weighted_Std(Current_X, :) = Weighted_Binned_Angle_Line_Weighted_Std;
    Weighted_Binned_Angle_Count(Current_X, :) = Weighted_Binned_Angle_Line_Count;
    disp(Current_X);
end
Parpool_Delete();
clear Weights X Y Weighted_Binned_Angle_Line_Min Weighted_Binned_Angle_Line_Max Weighted_Binned_Angle_Line_Mean Weighted_Binned_Angle_Line_Std Weighted_Binned_Angle_Line_Weighted_Mean Weighted_Binned_Angle_Line_Weighted_Std Weighted_Binned_Angle_Line_Count X_Index Y_Index;

figure();
%Surf_Fig = surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), Weighted_Binned_Angle_Min,'FaceAlpha', .8);
%set(Surf_Fig, 'linestyle', 'none');
imagesc(X_Bins, Y_Bins, Weighted_Binned_Angle_Min);
colorbar();
title('Min Angle');

figure();
%Surf_Fig = surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), Weighted_Binned_Angle_Max,'FaceAlpha', .8);
%set(Surf_Fig, 'linestyle', 'none');
imagesc(X_Bins, Y_Bins, Weighted_Binned_Angle_Max);
colorbar();
title('Max Angle');

figure();
%Surf_Fig = surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), Weighted_Binned_Angle_Mean,'FaceAlpha', .8);
%set(Surf_Fig, 'linestyle', 'none');
imagesc(X_Bins, Y_Bins, Weighted_Binned_Angle_Mean);
colorbar();
title('Mean Angle');

figure();
%Surf_Fig = surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), Weighted_Binned_Angle_Std,'FaceAlpha', .8);
%set(Surf_Fig, 'linestyle', 'none');
imagesc(X_Bins, Y_Bins, Weighted_Binned_Angle_Std);
colorbar();
title('Std Angle');

figure();
%Surf_Fig = surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), Weighted_Binned_Angle_Weighted_Mean,'FaceAlpha', .8);
%set(Surf_Fig, 'linestyle', 'none');
imagesc(X_Bins, Y_Bins, Weighted_Binned_Angle_Weighted_Mean);
colorbar();
title('Mean Weighted Angle');

figure();
%Surf_Fig = surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), Weighted_Binned_Angle_Weighted_Std,'FaceAlpha', .8);
%set(Surf_Fig, 'linestyle', 'none');
imagesc(X_Bins, Y_Bins, Weighted_Binned_Angle_Weighted_Std);
colorbar();
title('Std Weighted Angle');

figure();
%Surf_Fig = surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), Weighted_Binned_Angle_Count,'FaceAlpha', .8);
%set(Surf_Fig, 'linestyle', 'none');
imagesc(X_Bins, Y_Bins, Weighted_Binned_Angle_Count);
colorbar();
title('Count');

%% Weighted Histogram
Histogram_Bins = linspace(0, max(Event_Angle(:)) + 0.01, Num_Hist_Bins + 1);
Histogram_Angle_Min = zeros(1, length(Histogram_Bins) - 1);
Histogram_Angle_Max = zeros(1, length(Histogram_Bins) - 1);
Histogram_Angle_Mean = zeros(1, length(Histogram_Bins) - 1);
Histogram_Angle_Std = zeros(1, length(Histogram_Bins) - 1);
Histogram_Angle_Count = zeros(1, length(Histogram_Bins) - 1);
for Current_Histogram_Bin = 1:length(Histogram_Bins) - 1
    Index = (Histogram_Bins(Current_Histogram_Bin) <= Event_Angle) & (Event_Angle < Histogram_Bins(Current_Histogram_Bin + 1));
    if(sum(Index, 'omitnan') ~= 0)
        Angular_Weight = Pixel_Weight;
        Histogram_Angle_Min(1, Current_Histogram_Bin) = min(Angular_Weight);
        Histogram_Angle_Max(1, Current_Histogram_Bin) = max(Angular_Weight);
        Histogram_Angle_Mean(1, Current_Histogram_Bin) = mean(Angular_Weight, 'omitnan');
        Histogram_Angle_Std(1, Current_Histogram_Bin) = std(Angular_Weight, 'omitnan');
        Histogram_Angle_Count(1, Current_Histogram_Bin) = length(Angular_Weight);
    end
end

figure();
histogram('BinEdges', Histogram_Bins, 'BinCounts', Histogram_Angle_Min);
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Min Weight");

figure();
histogram('BinEdges', Histogram_Bins, 'BinCounts', Histogram_Angle_Max);
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Max Weight");

figure();
histogram('BinEdges', Histogram_Bins, 'BinCounts', Histogram_Angle_Mean);
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Mean Weight");

figure();
histogram('BinEdges', Histogram_Bins, 'BinCounts', Histogram_Angle_Std);
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Std Weight");

%Average histogram bins
Histogram_Bins_Average = (Histogram_Bins(1:end-1) + Histogram_Bins(2:end))./2;
figure();
errorbar(Histogram_Bins_Average, Histogram_Angle_Mean, Histogram_Angle_Std);

%% Weight-Energy histogram
Num_Hist_Bins = 150;
Num_Energy_Bins = 100;
Energy_Max = 100;
Histogram_Bins = linspace(0, max(Event_Angle(:)) + 0.01, Num_Hist_Bins + 1);
Energy_Bins = linspace(0, Energy_Max, Num_Energy_Bins + 1);
Energy_Histogram_Angle_Min = zeros(length(Energy_Bins) - 1, length(Histogram_Bins) - 1);
Energy_Histogram_Angle_Max = zeros(length(Energy_Bins) - 1, length(Histogram_Bins) - 1);
Energy_Histogram_Angle_Mean = zeros(length(Energy_Bins) - 1, length(Histogram_Bins) - 1);
Energy_Histogram_Angle_Std = zeros(length(Energy_Bins) - 1, length(Histogram_Bins) - 1);
Energy_Histogram_Angle_Count = zeros(length(Energy_Bins) - 1, length(Histogram_Bins) - 1);
for Current_Histogram_Bin = 1:length(Histogram_Bins) - 1
    Index = (Histogram_Bins(Current_Histogram_Bin) <= Event_Angle) & (Event_Angle < Histogram_Bins(Current_Histogram_Bin + 1));
    for Current_Energy_Bin = 1:length(Energy_Bins) - 1
        Energy_Index = (Energy_Bins(Current_Energy_Bin) <= Energy) & (Energy < Energy_Bins(Current_Energy_Bin + 1));
        Energy_Hist_Index = Index & Energy_Index;
        if(sum(Energy_Hist_Index, 'omitnan') ~= 0)
            Angular_Weight = Weight(Energy_Hist_Index);
            Energy_Histogram_Angle_Min(Current_Energy_Bin, Current_Histogram_Bin) = min(Angular_Weight);
            Energy_Histogram_Angle_Max(Current_Energy_Bin, Current_Histogram_Bin) = max(Angular_Weight);
            Energy_Histogram_Angle_Mean(Current_Energy_Bin, Current_Histogram_Bin) = mean(Angular_Weight, 'omitnan');
            Energy_Histogram_Angle_Std(Current_Energy_Bin, Current_Histogram_Bin) = std(Angular_Weight, 'omitnan');
            Energy_Histogram_Angle_Count(Current_Energy_Bin, Current_Histogram_Bin) = length(Angular_Weight);
        end
    end
end

figure();
imagesc(Histogram_Bins, Energy_Bins, Energy_Histogram_Angle_Min);
colorbar();
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Energy [keV]");
title("Min Weight");

figure();
imagesc(Histogram_Bins, Energy_Bins, Energy_Histogram_Angle_Max);
colorbar();
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Energy [keV]");
title("Max Weight");

figure();
imagesc(Histogram_Bins, Energy_Bins, Energy_Histogram_Angle_Mean);
colorbar();
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Energy [keV]");
title("Mean Weight");

figure();
imagesc(Histogram_Bins, Energy_Bins, Energy_Histogram_Angle_Std);
colorbar();
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Energy [keV]");
title("Std Weight");

figure();
imagesc(Histogram_Bins, Energy_Bins, Energy_Histogram_Angle_Count);
colorbar();
xlabel(strcat("Angular Incidence [", char(176), "]"));
ylabel("Energy [keV]");
title("Num Weight");
%%
%Create binned 2d histogram for x-y data
% xzCount = histcounts2(X(:), Y(:), X_Bins, Y_Bins);
% figure();
% surf(Grid_X, Grid_Y, min(zlim())*ones(size(Grid_X)), xzCount,'FaceAlpha', .8);



%Perform work on the MAT file
% for Current_Mat_File = 1:length(Mat_File_Path)
%     %Create new directory to place processed files into
%     [Output_Directory, Filename, Extension] = fileparts(Mat_File_Path{Current_Mat_File});
%     Parent_Directory = fileparts(Output_Directory);
%     Output_Directory = fullfile(Parent_Directory, 'Processed-Output');
%     Attempt_Directory_Creation(Output_Directory);
% 
%     %% Filter data within the MAT file
%     %Position the events land on the detection plane
% %     Filters.X.Min = -0.04;
% %     Filters.X.Max = 0.04;
% %     Filters.Y.Min = -0.04;
% %     Filters.Y.Max = 0.04;
% %     %Filters.Z.Min = 0;
% %     %Filters.Z.Max = 1;
% %     %Angle from the normal (to Z), must be +ve valued
% %     Filters.Angle.Min = 0;
% %     Filters.Angle.Max = 45;
% %     %Energy [KeV]
% %     Filters.Energy.Min = 0;
% %     Filters.Energy.Max = 130;
% %     %Weighting
% %     Filters.Weight.Min = 0.1;
%     %Filters.Weight.Max = 35;
%     Filters.Photons.Min = 0.05;
%     %Filters.Photons.Max = 35;
%     Filters.Photons.Interval = 0.01;
%     
%     %Create filepath for the filtered MAT file
%     Filtered_Mat_File_Path{Current_Mat_File} = fullfile(Output_Directory, strcat(Filename, '-Filtered', Extension));
%     %Filter the file according to parameters previously set
%     tic
%     Filtered_Mat_File = MCPL_Filter_MAT_Data(Mat_File_Path{Current_Mat_File}, Filtered_Mat_File_Path{Current_Mat_File}, Filters);
%     toc
%     
%     %Create filepath for the Recreated MCPL file
%     MCPL_Filepath = fullfile(Output_Directory, strcat(Filename, '-Processed.MCPL'));
%     %Convert MAT file back to an MCPL file
%     tic
%     MCPL_File = MAT_To_MCPL(Filtered_Mat_File, MCPL_Filepath, Display_Write_Progress);
%     toc
%     %Test extraction of MCPL file (only for comparison if not filtering data)
%     %Mat_File_Path_2 = MCPL_To_MAT(MCPL_File, Read_Parameters);
%     % Compare initial MAT and recreated MAT files match
%     %visdiff(Mat_File_Path{Current_Mat_File}, Mat_File_Path_2{1});
% end

% figure = Get_Figure([], true);
% for Current_Ray = 1:length(X)
%     figure = Get_Figure(figure, true);
%     hold on;
%     plot3([X(Current_Ray), X(Current_Ray)], [Y(Current_Ray), Y(Current_Ray)], [Z(Current_Ray), Z(Current_Ray)]);
% end
% figure = Get_Figure(figure, false);
