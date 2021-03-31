%% Filters a matfile containing MCPL data based on filter parameters specifed in Filters
function Filtered_Mat_File_Path = Filter_MPCL_MAT_Data(Mat_File_Path, Filtered_Mat_File_Path, Filters)
    %% Input handling
    if(nargin ~= 3)
        error("Filter_MCPL_MAT_Data: Expected 3 inputs.");
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
            error("Filter_MCPL_MAT_Data: Missing variables X, Y, Z or Weight within original MAT file");
        end
        if(all(ismember({'Header'}, Mat_File_Variables)))
%% TODO attempt to create header?
        else
            error("Filter_MCPL_MAT_Data: Missing MCPL Header within original MAT file");
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
                error("Filter_MCPL_MAT_Data: Error creating EKinDir in original MAT file.");
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
                error("Filter_MCPL_MAT_Data: Error creating Three Direction Vectors (Dx, Dy, Dz) and Energy in original MAT file.");
            end
        else
            error("Filter_MCPL_MAT_Data: Requires either [Three Direction Vectors (Dx, Dy, Dz) and Energy] or [EKinDir] for filtering.");
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
        if (~(std(Size_1(:)) == 0 && std(Size_2(:)) == 0))
            error("Filter_MCPL_MAT_Data: Mismatch in variable sizes.");
        end
        %If applicable; create the directory path to save the 
        Filtered_Mat_Directory_Path = fileparts(Filtered_Mat_File_Path);
        if(~isfolder(Filtered_Mat_Directory_Path))
            Directory_Created = Attempt_Directory_Creation(Filtered_Mat_Directory_Path);
            if(~Directory_Created)
                error(strcat("Filter_MCPL_MAT_Data: Could not create directory path to : ", Filtered_Mat_Directory_Path));
            end
        end
        %Check Filters is a structure
        if(~isstruct(Filters))
            error("Filter_MCPL_MAT_Data: Filters expected to be a structure.");
        end
        %Verify the fields in the structure against valid filter types
        Allowed_Filters = {'X', 'Y', 'Z', 'Energy', 'Weight', 'Angle'};
        %Optional filters (based on data available in the original MAT file)
        if(Exists.Polarisation)
            Allowed_Filters(end + 1: end + 3) = {'Px', 'Py' ,'Px'};
        else
            disp("Filter_MCPL_MAT_Data: Polarisation filtering disabled, no relevant data in MAT file");
        end
        if(Exists.Userflag)
            Allowed_Filters(end + 1) = {'UserFlag'};
            disp("Filter_MCPL_MAT_Data: Userflag filtering disabled, no relevant data in MAT file");
        end
        if(Exists.Time)
            Allowed_Filters(end + 1) = {'Time'};
            disp("Filter_MCPL_MAT_Data: Time filtering disabled, no relevant data in MAT file");
        end
        Active_Filters = fieldnames(Filters);
        Remove_Filters = zeros(size(Active_Filters), 'logical');
        for Current_Active_Filter = 1:length(Active_Filters)
            Filter_Allowed = any(strcmpi(Allowed_Filters, Active_Filters{Current_Active_Filter}));
            if(~Filter_Allowed)
                Remove_Filters(Current_Active_Filter) = 1;
                disp(strcat("Filter_MCPL_MAT_Data: Ignoring Invalid Filter : ", Active_Filters{Current_Active_Filter}));
            end
        end
        %Delete any invalid filter requests
        Active_Filters(Remove_Filters) = [];
        
        %% Create a list of all indexed present within the original MAT file
        Allowed_Index_List = ones(Size_1(1), 1, 'logical');
        
        %% Generic sense-checks for structure inputs, then retain indicies of valid data meeting all criteria
        for Current_Active_Filter = 1:length(Active_Filters)
            if(isstruct(Filters.(Active_Filters{Current_Active_Filter})))
                %Find Min and Max fields
                Filter_Min_Max = fieldnames(Filters.(Active_Filters{Current_Active_Filter}));
                Filters.(Active_Filters{Current_Active_Filter}).Min_Active = false;
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
                %Minimum value exceptions for angle end energy
                if(strcmpi(Active_Filters{Current_Active_Filter}, 'Energy') || strcmpi(Active_Filters{Current_Active_Filter}, 'Angle'))
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
                            warning(strcat("Filter_MCPL_MAT_Data : Minimum ", Active_Filters{Current_Active_Filter}, " trunicated to 90"));
                        end
                    end
                end
                % load data and set Allowed_Index_List to false where the current event doesn't meet the required conditions.
                if(Filters.(Active_Filters{Current_Active_Filter}).Max_Active || Filters.(Active_Filters{Current_Active_Filter}).Min_Active)
                    if(strcmpi(Active_Filters{Current_Active_Filter}, 'Angle'))
                        %load Calculate angle from Z dimension (0,0,1) due to Dx, Dy, Dz
                        Data = acosd((Mat_File_Reference.Dz)./sqrt(Mat_File_Reference.Dx.^2 + Mat_File_Reference.Dy.^2 + Mat_File_Reference.Dz.^2));
                    else
                        %Load data specific to the relative variable and cut based on min/max thresholds
                        Data = Mat_File_Reference.(Active_Filters{Current_Active_Filter});
                    end
                    %Compare data to limits, set indexes to false when out of range (min / max)
                    if(Filters.(Active_Filters{Current_Active_Filter}).Max_Active)
                        Remove_List = Data > Filters.(Active_Filters{Current_Active_Filter}).Max;
                        disp(strcat("Filter_MCPL_MAT_Data : Removed ", num2str(sum(Remove_List)), " Events due to Max ", Active_Filters{Current_Active_Filter})); 
                        Allowed_Index_List(Remove_List) = false;
                    end
                    if(Filters.(Active_Filters{Current_Active_Filter}).Min_Active)
                        Remove_List = Data < Filters.(Active_Filters{Current_Active_Filter}).Min;
                        Allowed_Index_List(Remove_List) = false;
                        disp(strcat("Filter_MCPL_MAT_Data : Removed ", num2str(sum(Remove_List)), " Events due to Min ", Active_Filters{Current_Active_Filter}));
                    end
                end
            else
                disp(strcat("Filter_MCPL_MAT_Data: Ignoring Invalid Filter Structure : ", Active_Filters{Current_Active_Filter}));
            end
        end
        clear Filter_Min_Max Current_Active_Filter;
        clear Size_1 Size_2 Filtered_Mat_Directory_Path Remove_Filters Data;
        
        %Create filtered MAT file
        Filtered_Mat_File_Reference = matfile(Filtered_Mat_File_Path);
        Filtered_Mat_File_Reference.Properties.Writable = true;
        %Copy valid data
        for Current_File_Row = 1:length(Allowed_Index_List)
            
        end
    else
        error("Filter_MCPL_Mat_Data: MAT file not found");
    end
end