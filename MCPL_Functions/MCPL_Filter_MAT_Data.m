%% Filters a matfile containing MCPL data based on filter parameters specifed in Filters
function Filtered_Mat_File_Path = MCPL_Filter_MAT_Data(Mat_File_Path, Filtered_Mat_File_Path, Filters)
    %% Input handling
    if(nargin ~= 3)
        error("Filter_MCPL_MAT_Data : Expected 3 inputs.");
    end
    if(isfile(Mat_File_Path))
        %Open reference to the MAT file containing the unfiltered data
        Mat_File_Reference = matfile(Mat_File_Path);
        %Identify variables within the MAT file
        try
            Mat_File_Variables = who('-file', Mat_File_Path);
        catch
            error(strcat("Invalid MAT file format: ", Mat_File_Path));
        end
        %% Validate relevant variables exist within the MAT fil
        %X, Y, Z, Weight
        if(all(ismember({'X','Y','Z','Weight'}, Mat_File_Variables)))
        else
            error("Filter_MCPL_MAT_Data : Missing variables X, Y, Z or Weight within original MAT file");
        end
        if(all(ismember({'Header'}, Mat_File_Variables)))
            %Load Header
            Header = Mat_File_Reference.Header;
        else
            error("Filter_MCPL_MAT_Data : Missing MCPL Header within original MAT file");
        end
        %Dx, Dy, Dz
        if(all(ismember({'Dx', 'Dy', 'Dz'}, Mat_File_Variables)))
            Exists.Vectors = true;
        else
            Exists.Vectors = false;
        end
        %Energy
        if(all(ismember({'Energy'}, Mat_File_Variables)))
            Exists.Energy = true;
        else
            Exists.Energy = false;
        end
        %EKinDir
        if(all(ismember({'EKinDir_1','EKinDir_2','EKinDir_3'}, Mat_File_Variables)))
            Exists.EKinDir = true;
        else
            Exists.EKinDir = false;
        end
        %Validate EKinDir / Energy + Direction Vectors exist, recreate
        if(Exists.Energy && Exists.Vectors && Exists.EKinDir)
            %Nothing required
        elseif(Exists.Energy && Exists.Vectors && ~Exists.EKinDir)
            %Recalculate EKinDir
            Mat_File_Reference.Properties.Writable = true;
            [Mat_File_Reference.EKinDir_1, Mat_File_Reference.EKinDir_2, Mat_File_Reference.EKinDir_3] = EKinDir_Pack(Mat_File_Reference.Dx, Mat_File_Reference.Dy, Mat_File_Reference.Dz, Mat_File_Reference.Energy);
            Mat_File_Reference.Properties.Writable = false;
            %Update list of variables
            Mat_File_Variables = who('-file', Mat_File_Path);
            %Re-check existance of newly added variables
            %EKinDir
            if(all(ismember({'EKinDir_1','EKinDir_2','EKinDir_3'}, Mat_File_Variables)))
                Exists.EKinDir = true;
            else
                Exists.EKinDir = false;
            end
            if(~Exists.EKinDir)
                error("Filter_MCPL_MAT_Data : Error creating EKinDir in original MAT file.");
            end
        elseif(~Exists.Energy && ~Exists.Vectors && Exists.EKinDir)
            %Recalculate Energy and Direction Vectors
            Mat_File_Reference.Properties.Writable = true;
            [Mat_File_Reference.Dx, Mat_File_Reference.Dy, Mat_File_Reference.Dz, Mat_File_Reference.Energy] = EKinDir_Unpack(Mat_File_Reference.EKinDir_1, Mat_File_Reference.EKinDir_2, Mat_File_Reference.EKinDir_3, Mat_File_Reference.Header.MCPL_Version);
            Mat_File_Reference.Properties.Writable = false;
            %Update list of variables
            Mat_File_Variables = who('-file', Mat_File_Path);
            %Re-check existance of newly added variables
            %Energy
            if(all(ismember({'Energy'}, Mat_File_Variables)))
                Exists.Energy = true;
            else
                Exists.Energy = false;
            end
            %Dx, Dy, Dz
            if(all(ismember({'Dx', 'Dy', 'Dz'}, Mat_File_Variables)))
                Exists.Vectors = true;
            else
                Exists.Vectors = false;
            end
            if(~Exists.EKinDir || ~Exists.Vectors)
                error("Filter_MCPL_MAT_Data : Error creating Three Direction Vectors (Dx, Dy, Dz) and Energy in original MAT file.");
            end
        else
            error("Filter_MCPL_MAT_Data : Requires either [Three Direction Vectors (Dx, Dy, Dz) and Energy] or [EKinDir] for filtering.");
        end
        if(all(ismember({'Px','Py','Pz'}, Mat_File_Variables)))
            Exists.Polarisation = true;
        else
            Exists.Polarisation = false;
        end
        if(all(ismember({'UserFlag'}, Mat_File_Variables)))
            Exists.Userflag = true;
        else
            Exists.Userflag = false;
        end
        if(all(ismember({'Time'}, Mat_File_Variables)))
            Exists.Time = true;
        else
            Exists.Time = false;
        end
        %Remove header from the variable list (it's only ever a 1x1 struct)
        Mat_File_Variables(strcmp(Mat_File_Variables, 'Header')) = [];
        %Validate the length of each field in the structure matches
        Size_1 = zeros(length(Mat_File_Variables),1);
        Size_2 = zeros(length(Mat_File_Variables),1);
        for Current_Field = 1:length(Mat_File_Variables)
            [Size_1(Current_Field), Size_2(Current_Field)] = size(Mat_File_Reference, Mat_File_Variables{Current_Field});
        end
        if (~(range(Size_1(:)) == 0 && range(Size_2(:)) == 0))
            error("Filter_MCPL_MAT_Data : Mismatch in variable sizes.");
        end
        % Verify the length of the fields matches the number of events in the header
        if(~(Floating_Point_Equal(mean(Size_1(:)), Header.Particles) || Floating_Point_Equal(mean(Size_2(:)), Header.Particles)))
            error("Filter_MCPL_MAT_Data : Mismatch in variable sizes compared to number of total events in the header.");
        end
        %If applicable; create the directory path to save the 
        Filtered_Mat_Directory_Path = fileparts(Filtered_Mat_File_Path);
        if(~isfolder(Filtered_Mat_Directory_Path))
            Directory_Created = Attempt_Directory_Creation(Filtered_Mat_Directory_Path);
            if(~Directory_Created)
                error(strcat("Filter_MCPL_MAT_Data : Could not create directory path to : ", Filtered_Mat_Directory_Path));
            end
        end
        %Check Filters is a structure
        if(~isstruct(Filters))
            error("Filter_MCPL_MAT_Data: Filters expected to be a structure.");
        end
        %Verify the fields in the structure against valid filter types
        Allowed_Filters = {'X', 'Y', 'Z', 'Energy', 'Weight', 'Photons', 'Angle'};
        %Optional filters (based on data available in the original MAT file)
        if(Exists.Polarisation)
            Allowed_Filters(end + 1: end + 3) = {'Px', 'Py' ,'Px'};
        else
            disp("Filter_MCPL_MAT_Data : Polarisation filtering disabled, no relevant data in MAT file");
        end
        if(Exists.Userflag)
            Allowed_Filters(end + 1) = {'UserFlag'};
            disp("Filter_MCPL_MAT_Data : Userflag filtering disabled, no relevant data in MAT file");
        end
        if(Exists.Time)
            Allowed_Filters(end + 1) = {'Time'};
            disp("Filter_MCPL_MAT_Data : Time filtering disabled, no relevant data in MAT file");
        end
        Active_Filters = fieldnames(Filters);
        Remove_Filters = zeros(size(Active_Filters), 'logical');
        for Current_Active_Filter = 1:length(Active_Filters)
            Filter_Allowed = any(strcmpi(Allowed_Filters, Active_Filters{Current_Active_Filter}));
            if(~Filter_Allowed)
                Remove_Filters(Current_Active_Filter) = 1;
                disp(strcat("Filter_MCPL_MAT_Data : Ignoring Invalid Filter : ", Active_Filters{Current_Active_Filter}));
            end
        end
        %Delete any invalid filter requests
        Active_Filters(Remove_Filters) = [];
        
        %% Generic sense-checks for structure inputs, then retain indicies of valid data meeting all criteria
        disp("Filter_MCPL_MAT_Data : Validating Filtering Inputs.");
        %% Load datastore file (excluding the header)
        File_Data_Store = tall(fileDatastore(Mat_File_Path, 'ReadFcn', @(x)struct2table(load(x,'-regexp', '^(?!Header$).')), 'UniformRead', true));
        File_Variable_Names = gather(File_Data_Store.Properties.VariableNames);
        disp("Filter_MCPL_MAT_Data : Finding number of events pre-filtering.")
        Num_Events_Start = gather(height(File_Data_Store));
        disp("Filter_MCPL_MAT_Data : Note one event can be filtered out by multiple individual conditions.");
        for Current_Active_Filter = 1:length(Active_Filters)
            if(isstruct(Filters.(Active_Filters{Current_Active_Filter})))
                %% Validate filter input
                %Find Min and Max fields for current filter
                Filter_Min_Max = fieldnames(Filters.(Active_Filters{Current_Active_Filter}));
                %Verify Min field exists
                if(any(strcmpi(Filter_Min_Max, 'Min')))
                    if(isnumeric(Filters.(Active_Filters{Current_Active_Filter}).Min))
                        Filters.(Active_Filters{Current_Active_Filter}).Min_Active = true;
                    else
                        warning(strcat("Filter_MCPL_MAT_Data : Expected Numeric Input for Minimum ", Active_Filters{Current_Active_Filter}, " : Ignoring Filter"));
                    end
                end
                Filters.(Active_Filters{Current_Active_Filter}).Max_Active = false;
                %Verify Max field exists
                if(any(strcmpi(Filter_Min_Max, 'Max')))
                    if(isnumeric(Filters.(Active_Filters{Current_Active_Filter}).Max))
                        Filters.(Active_Filters{Current_Active_Filter}).Max_Active = true;
                    else
                        warning(strcat("Filter_MCPL_MAT_Data : Expected Numeric Input for Maximum ", Active_Filters{Current_Active_Filter}, " : Ignoring Filter"));
                    end
                end
                %If both min and max are active; check they are appropriately ordered (min < max)
                if(Filters.(Active_Filters{Current_Active_Filter}).Max_Active && Filters.(Active_Filters{Current_Active_Filter}).Min_Active)
                    if(Filters.(Active_Filters{Current_Active_Filter}).Max < Filters.(Active_Filters{Current_Active_Filter}).Min)
                        Temp = Filters.(Active_Filters{Current_Active_Filter}).Max;
                        Filters.(Active_Filters{Current_Active_Filter}).Max = Filters.(Active_Filters{Current_Active_Filter}).Min;
                        Filters.(Active_Filters{Current_Active_Filter}).Min = Temp;
                        warning(strcat("Filter_MCPL_MAT_Data : Swapped Min and Max bounds for ", Active_Filters{Current_Active_Filter}));
                    end
                end
                %Minimum value exceptions for angle, energy, photons
                if(any(strcmpi(Active_Filters{Current_Active_Filter}, {'Energy', 'Angle', 'Photons'})))
                    if(Filters.(Active_Filters{Current_Active_Filter}).Min_Active)
                        if(Filters.(Active_Filters{Current_Active_Filter}).Min < 0)
                            Filters.(Active_Filters{Current_Active_Filter}).Min = 0;
                            warning(strcat("Filter_MCPL_MAT_Data : Minimum ", Active_Filters{Current_Active_Filter}, " trunicated to 0"));
                        end
                    end
                end
                %maximum value exception for angle
                if(strcmpi(Active_Filters{Current_Active_Filter}, 'Angle'))
                    if(Filters.(Active_Filters{Current_Active_Filter}).Max_Active)
                        if(Filters.(Active_Filters{Current_Active_Filter}).Max > 90)
                            Filters.(Active_Filters{Current_Active_Filter}).Max = 90;
                            warning(strcat("Filter_MCPL_MAT_Data : Maximum ", Active_Filters{Current_Active_Filter}, " trunicated to 90"));
                        end
                    end
                end
                %% Filter Data using current filter (data that directly exists in the table)
                %Verify corresponding field exists directly in data
                if(any(strcmpi(Active_Filters(Current_Active_Filter), File_Variable_Names)))
                    %Minimum filter
                    if(Filters.(Active_Filters{Current_Active_Filter}).Min_Active)
                        Remove_Index = File_Data_Store.(Active_Filters{Current_Active_Filter}) < Filters.(Active_Filters{Current_Active_Filter}).Min;
                        File_Data_Store(Remove_Index,:) = [];
                        disp(strcat("Filter_MCPL_MAT_Data : Removed ", num2str(gather(sum(Remove_Index))), " Events due to Min ", Active_Filters{Current_Active_Filter}));
                    end
                    %Maximum filter
                    if(Filters.(Active_Filters{Current_Active_Filter}).Max_Active)
                        Remove_Index = File_Data_Store.(Active_Filters{Current_Active_Filter}) > Filters.(Active_Filters{Current_Active_Filter}).Max;
                        File_Data_Store(Remove_Index,:) = [];
                        disp(strcat("Filter_MCPL_MAT_Data : Removed ", num2str(gather(sum(Remove_Index))), " Events due to Max ", Active_Filters{Current_Active_Filter}));
                    end
                else
                    %% Data that requires prior calculation from original table data to filter
                    if(any(strcmpi(Active_Filters{Current_Active_Filter}, {'Angle', 'Photons'})))
                        %% Angular deviation from Z directional vector
                        if(strcmpi(Active_Filters{Current_Active_Filter}, 'Angle'))
                            %Calculate angle from Z dimension (0,0,1) due to Dx, Dy, Dz
                            Calculated_Data = acosd((File_Data_Store.Dz)./sqrt(File_Data_Store.Dx.^2 + File_Data_Store.Dy.^2 + File_Data_Store.Dz.^2));
                            %Minimum filter
                            if(Filters.(Active_Filters{Current_Active_Filter}).Min_Active)
                                Remove_Index = Calculated_Data < Filters.(Active_Filters{Current_Active_Filter}).Min;
                                File_Data_Store(Remove_Index,:) = [];
                                disp(strcat("Filter_MCPL_MAT_Data : Removed ", num2str(gather(sum(Remove_Index))), " Events due to Min ", Active_Filters{Current_Active_Filter}));
                            end
                            %Maximum filter
                            if(Filters.(Active_Filters{Current_Active_Filter}).Max_Active)
                                Remove_Index = Calculated_Data > Filters.(Active_Filters{Current_Active_Filter}).Max;
                                File_Data_Store(Remove_Index,:) = [];
                                disp(strcat("Filter_MCPL_MAT_Data : Removed ", num2str(gather(sum(Remove_Index))), " Events due to Max ", Active_Filters{Current_Active_Filter}));
                            end
                        end
                        %% By photon contribution at specific p-value bands
                        if(strcmpi(Active_Filters{Current_Active_Filter}, 'Photons'))
                            %% Get contribution of each group of p-values in terms of photons
                            %Get histogram data (hidden figure)
                            Histogram_Figure = Get_Figure([], true);
                            Histogram_Num_Bins = round(sqrt(gather(height(File_Data_Store))));
                            Histogram = histogram(File_Data_Store.Weight, Histogram_Num_Bins);
                            %Get relevant data from the histogram figure
                            Histogram_Bin_Edges = Histogram.BinEdges;
                            Histogram_Bin_Width = Histogram.BinWidth;
                            Histogram_Counts = Histogram.Values;
                            %Close the hidden histogram figure
                            close(Histogram_Figure);
                            %get intervals for the histogram bins in a 1:1 index relationship
                            Histogram_Bin_Start = Histogram_Bin_Edges(1:end-1);
                            Histogram_Bin_End = Histogram_Bin_Edges(2:end);
                            Histogram_Bins = Histogram_Bin_Start + Histogram_Bin_Width/2;
                            Photon_Contribution = Histogram_Bins.*Histogram_Counts;
                            Relative_Photon_Contribution = Photon_Contribution./max(Photon_Contribution);
                            Absolute_Photon_Contribution = Photon_Contribution./sum(Photon_Contribution);
                            %% find groups of sequential histogram bins to remove
                            Remove_Histogram_Bins = zeros(size(Relative_Photon_Contribution),'logical');
                            if(Filters.(Active_Filters{Current_Active_Filter}).Min_Active)
                                Remove_Histogram_Bins = Remove_Histogram_Bins | (Relative_Photon_Contribution < Filters.(Active_Filters{Current_Active_Filter}).Min);
                            end
                            if(Filters.(Active_Filters{Current_Active_Filter}).Max_Active)
                                Remove_Histogram_Bins = Remove_Histogram_Bins | (Relative_Photon_Contribution > Filters.(Active_Filters{Current_Active_Filter}).Max);
                            end
                            [Group_Index_Start, Group_Index_End] = Find_Logical_Groups(Remove_Histogram_Bins);
                            if(numel(Group_Index_Start) ~= numel(Group_Index_End))
                                error("Filter_MCPL_MAT_Data : Expected photon groups to have the same number of starting and ending indicies.");
                            end
                            
                            %% Filtering
                            %Create a 1:1 index list corresponding to the histogram
                            Remove_Index = File_Data_Store.X;
                            Remove_Index(:) = 0;
                            %Filter out by group
                            Total_Percent_Filtered = 0;
                            for Current_Group_Index = 1:length(Group_Index_Start)
                                Total_Percent_Filtered = Total_Percent_Filtered + sum(Absolute_Photon_Contribution(Group_Index_Start:Group_Index_End))*100;
                                disp(strcat("Filter_MCPL_MAT_Data : Filtering photon group ",num2str(Current_Group_Index), " between the p interval ", num2str(Histogram_Bin_Start(Current_Group_Index))," to ", num2str(Histogram_Bin_End(Current_Group_Index)), " with a contribution of ", num2str(sum(Absolute_Photon_Contribution(Group_Index_Start(Current_Group_Index):Group_Index_End(Current_Group_Index)))*100), "% of photons."));
                                Remove_Index = Remove_Index | ((Histogram_Bin_Start(Current_Group_Index) >= File_Data_Store.Weight) & (File_Data_Store.Weight < Histogram_Bin_End(Current_Group_Index)));
                            end
                            File_Data_Store(Remove_Index,:) = [];
                            disp(strcat("Filter_MCPL_MAT_Data : Removed ", num2str(gather(sum(Remove_Index))), " Events due to ", Active_Filters{Current_Active_Filter}));
                        end
                    else
                        disp(strcat("Filter_MCPL_MAT_Data : No corresponding data found for filtering: ", Active_Filters{Current_Active_Filter}));
                    end
                end
            else
                disp(strcat("Filter_MCPL_MAT_Data: Ignoring Invalid Filter Structure: ", Active_Filters{Current_Active_Filter}));
            end
        end
        
        %% Valid data remains; display final output
        %Get the number of data remaining
        Num_Events_End = gather(height(File_Data_Store));

        %% Output
        %Check valid data remains
        if(Num_Events_End == 0)
            warning("Filter_MCPL_MAT_Data : No data remains after filtering, skipping writing file. Returning raw data file as output.");
            Filtered_Mat_File_Path = Mat_File_Path;
        else
            %Write processed datastore data to files
            Datastore_Directory_Path = fullfile(fileparts(Filtered_Mat_File_Path), 'Filtered');
            %Ensure the output directory is empty
            if(isfolder(Datastore_Directory_Path))
                disp("MCPL_To_MAT : Datastore directory already exists, clearing directory contents.");
                Attempt_Directory_Deletion(Datastore_Directory_Path);
                Attempt_Directory_Creation(Datastore_Directory_Path);
            end
            disp("Filter_MCPL_MAT_Data : Performing Datastore Operations and Saving to Datastore Partitions.");
            write(fullfile(Datastore_Directory_Path, 'Partition_*.mat'), File_Data_Store, 'WriteFcn', @Write_Data);
            %% Merge datastore files
            Filtered_Mat_File_Path = MCPL_Merge_Chunks(Datastore_Directory_Path, Header, Filtered_Mat_File_Path, true);
            %% Display output file progress
            disp(strcat("Filter_MCPL_MAT_Data : Input Events    : ", num2str(Num_Events_Start)));
            disp(strcat("Filter_MCPL_MAT_Data : Removed Events  : ", num2str(Num_Events_Start - Num_Events_End)));
        end
    else
        error("Filter_MCPL_Mat_Data : MAT file not found.");
    end
    
    %% Datastore write function (local to parent to save on memory duplication)
    function Write_Data(info, data)
        %Turn table into structure
        data = table2struct(data);
        %Turn vector structure into a scalar structure containing vector data
        data = Structure_Vector_To_Scalar(data);
        data.Datastore_Partition_Information = info;
        %Write data to file
        save(info.SuggestedFilename, '-v7.3', '-struct', 'data');
    end
end