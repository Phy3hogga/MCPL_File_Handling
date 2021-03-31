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
        %% Validate relevant variables exist within the MAT file
        if(all(ismember({'X','Y','Z','Dx', 'Dy', 'Dz','Weight'}, Mat_File_Variables)))
        else
            error("Filter_MCPL_MAT_Data: Missing key variables within original MAT file");
        end
        if(all(ismember({'Header'}, Mat_File_Variables)))
        else
            error("Filter_MCPL_MAT_Data: Missing key variables within original MAT file");
        end
        if(all(ismember({'EKinDir_1','EKinDir_2','EKinDir_3'}, Mat_File_Variables)))
            Exists.EKinDir= true;
        else
            Exists.EKinDir = false;
        end
        if(all(ismember({'Energy'}, Mat_File_Variables)))
            Exists.Energy = true;
        else
            Exists.Energy = false;
        end
        %Require either EKinDir or Energy to exist
        if((Exists.EKinDir == false) && (Exists.Energy == false))
            error("Filter_MCPL_MAT_Data: Missing Both EKinDir and Energy, one required.");
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
                error(strcat("Could not create directory path to : ", Filtered_Mat_Directory_Path));
            end
        end
        %Check Filters is a structure
        if(~isstruct(Filters))
            error("Filter_MCPL_MAT_Data: Filters expected to be a structure.");
        end
        %Verify the fields in the structure against valid filter types
        Allowed_Filters = {'X', 'Y', 'Z', 'Energy', 'Weight', 'Angle'};
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
        
        %% Generic sense-checks for structure inputs
        %Default values
        for Current_Active_Filter = 1:length(Active_Filters)
            Filters.(Active_Filters{Current_Active_Filter}).Min_Active = false;
            Filters.(Active_Filters{Current_Active_Filter}).Max_Active = false;
        end
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
                    Data = Mat_File_Reference.({Current_Active_Filter});
                    
                end
            else
                disp(strcat("Filter_MCPL_MAT_Data: Ignoring Invalid Filter Structure : ", Active_Filters{Current_Active_Filter}));
            end
        end
        clear Filter_Min_Max Current_Active_Filter;
        
        
        %% Create a list of all indexed present within the original MAT file
        Allowed_Index_List = ones(Size_1(1), 1, 'logical');
        clear Size_1 Size_2 Filtered_Mat_Directory_Path Remove_Filters;
        
        %Create filtered MAT file
        Filtered_Mat_File_Reference = matfile(Filtered_Mat_File_Path);
        Filtered_Mat_File_Reference.Properties.Writable = true;
    else
        error("Filter_MCPL_Mat_Data: MAT file not found");
    end
end