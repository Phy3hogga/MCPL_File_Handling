%% Merge MCPL File Chunks into a single MAT file
function Merged_File_Path = MCPL_Merge_Files(Datastore_Directory_Path, Output_File_Path, Remove_Temp_Files)
    %File path to write the merged data to
    [Root_File_Path, Filename, ~] = fileparts(Output_File_Path);
    Merged_File_Path = fullfile(Root_File_Path, strcat(Filename, ".mat"));
    
    MAT_File_List = Search_Files(Datastore_Directory_Path, '.mat');
    
    %% Recombination of the different partitions
    disp("MCPL_Merge_Files : Merging Datastore Partitions.");
    %Load matfile references to each file containing a chunk of processed data
    Chunk_File_Paths = fullfile({MAT_File_List(:).folder}, {MAT_File_List(:).name});
    Chunk_Matfile_References = cellfun(@matfile, Chunk_File_Paths, 'UniformOutput', false);
    %Find limits of array sizes for each chunk
    for Current_File = 1:length(Chunk_Matfile_References)
        Header = Chunk_Matfile_References{Current_File}.Header;
        Chunk_Events(Current_File) = Header.Particles;
    end
    Total_Chunk_Events = sum(Chunk_Events(:));

    %Create merged MAT file
    Merged_File_Reference = matfile(Merged_File_Path);
    Merged_File_Reference.Properties.Writable = true;
    
    %Preallocate variables within the merged MAT file
    if(Header.Opt_SinglePrecision)
        Empty_Byte_Type = single(0);
    else
        Empty_Byte_Type = double(0);
    end
    if(Header.Opt_Polarisation)
        Merged_File_Reference.Px(Total_Chunk_Events, 1) = Empty_Byte_Type;
        Merged_File_Reference.Py(Total_Chunk_Events, 1) = Empty_Byte_Type;
        Merged_File_Reference.Pz(Total_Chunk_Events, 1) = Empty_Byte_Type;
    end
    Merged_File_Reference.X(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.Y(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.Z(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.Dx(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.Dy(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.Dz(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.Energy(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.Time(Total_Chunk_Events, 1) = Empty_Byte_Type;
    if(Header.Save_EKinDir)
        Merged_File_Reference.EKinDir_1(Total_Chunk_Events, 1) = Empty_Byte_Type;
        Merged_File_Reference.EKinDir_2(Total_Chunk_Events, 1) = Empty_Byte_Type;
        Merged_File_Reference.EKinDir_3(Total_Chunk_Events, 1) = Empty_Byte_Type;
    end
    if(~Header.Opt_UniversalWeight)
        Merged_File_Reference.Weight(Total_Chunk_Events, 1) = Empty_Byte_Type;
    end
    if(Header.Opt_UniversalPDGCode == 0)
        Merged_File_Reference.PDGCode(Total_Chunk_Events, 1) = int32(0);
    end
    if(Header.Opt_Userflag)
        Merged_File_Reference.UserFlag(Total_Chunk_Events, 1) = uint32(0);
    end
    clear Empty_Byte_Type;

    %% Combine the different partitions / raw file data
    Remove_Variables = {'Datastore_Partition_Information', 'Header'};
    %Index of current row within the file to write to
    File_Write_Index = 1;
    for Current_File = 1:length(Chunk_Matfile_References)
        %Get variables in file
        Mat_File_Variables = whos(Chunk_Matfile_References{Current_File});
        %Remove any unwanted variables
        if(~isempty(Remove_Variables))
            for Current_Variable = 1:length(Remove_Variables)
                Mat_File_Variables(strcmpi({Mat_File_Variables.name}, Remove_Variables{Current_Variable})) = [];
            end
        end
        %Final index to write to
        File_Write_Index_End = File_Write_Index + Chunk_Events(Current_File) - 1;
        %Copy data for each variable to the combined file
        for Current_Variable = 1:length(Mat_File_Variables)
            %Get size of current variable
            Variable_Size = Mat_File_Variables(Current_Variable).size;
            if(~isempty(Variable_Size))
                %Copy data
                Merged_File_Reference.(Mat_File_Variables(Current_Variable).name)(File_Write_Index:File_Write_Index_End, 1) = reshape(Chunk_Matfile_References{Current_File}.(Mat_File_Variables(Current_Variable).name)(1:Variable_Size(1), 1:Variable_Size(2)), [], 1);
            else
                warning(strcat("MCPL_Merge_Chunks : Could not append data for variable : ", Mat_File_Variables(Current_Variable), " from partition ", Chunk_File_Paths(Current_Variable)));
            end
        end
        %Increment for next write pass
        File_Write_Index = File_Write_Index_End + 1;
    end
    
    %% Edit header information with removed events
    %Edit the number of events in the stored file
    Header.Particles = File_Write_Index_End;
    %Remove File_Chunks from the header data (No longer needed now data is combined)
    Header.File_Chunks = [];
    
    %% Copy header into the data file last (ensures writing is finished; if missing, file is invalid)
    Merged_File_Reference.Header = Header;
    
    %% Delete datastore directory (including contents)
    if(Remove_Temp_Files)
        Attempt_Directory_Deletion(Datastore_Directory_Path);
    end
end