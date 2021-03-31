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
        %Validate relevant variables exist within the MAT file
        if(all(ismember({'Header','X','Y','Z','Dx', 'Dy', 'Dz','Weight'}, Mat_File_Variables)))
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
        %If applicable; create the directory path to save the 
        Filtered_Mat_File_Path = fileparts(Filtered_Mat_File_Path);
        if(~isfolder(Filtered_Mat_File_Path))
            Directory_Created = Attempt_Directory_Creation(Filtered_Mat_File_Path);
            if(~Directory_Created)
                error(strcat("Could not create directory tree to : ", Filtered_Mat_File_Path));
            end
        end
        %Filtered_Mat_File_Reference = matfile();
    else
        error("Filter_MCPL_Mat_Data: MAT file not found");
    end
end