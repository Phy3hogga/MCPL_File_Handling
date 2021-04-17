%% Order datastore partition files
function File_Chunks = Order_Datastore_Partitions(Datastore_Directory_Path)
    %% Validate and Order the datastore partition files
    %Find partition files in the Datastore directory
    Files = Search_Files(Datastore_Directory_Path, '.mat');
    %Track number of partitions, partition index and variable length for each file
    Number_Of_Partitions = zeros(size(Files));
    Partition_Index = zeros(size(Files));
    Data_Length = zeros(size(Files));
    %Open references to all mat files
    Matfile_Partition_References = cellfun(@matfile, fullfile({Files.folder}, {Files.name}), 'UniformOutput', false);
    for Current_Reference = 1:length(Matfile_Partition_References)
        %Check datastore partition information exists in the partition file, load
        Mat_File_Variables = whos(Matfile_Partition_References{Current_Reference});
        if(~contains('Datastore_Partition_Information', {Mat_File_Variables.name}))
            error("MCPL_To_MAT : Datastore_Partition_Information missing from Datastore Partition");
        end
        Datastore_Partition_Information = Matfile_Partition_References{Current_Reference}.Datastore_Partition_Information;
        %Remove datastore from the variable list in memory as it's now been verfied to exist and has been read
        Mat_File_Variables(strcmp({Mat_File_Variables.name}, 'Datastore_Partition_Information')) = [];
        %Get number of total partitions
        if(~isfield(Datastore_Partition_Information, 'NumberOfPartitions'))
            error("MCPL_To_MAT : Datastore_Partition_Information missing NumberOfPartitions from Datastore Partition Information");
        end
        Number_Of_Partitions(Current_Reference) = Datastore_Partition_Information.NumberOfPartitions;
        %Get current partition
        if(~isfield(Datastore_Partition_Information, 'PartitionIndex'))
            error("MCPL_To_MAT : Datastore_Partition_Information missing PartitionIndex from Datastore Partition Information");
        end
        Partition_Index(Current_Reference) = Datastore_Partition_Information.PartitionIndex;
        %Get length of data in all fields
        Size_1 = zeros(length(Mat_File_Variables),1);
        Size_2 = zeros(length(Mat_File_Variables),1);
        for Current_Field = 1:length(Mat_File_Variables)
            Size_1(Current_Field) = Mat_File_Variables(Current_Field).size(1);
            Size_2(Current_Field) = Mat_File_Variables(Current_Field).size(2);
        end
        if (~(range(Size_1(:)) == 0 && range(Size_2(:)) == 0))
            error("MCPL_To_MAT : Mismatch in variable sizes.");
        end
        if(mean(Size_1(:)) == 1)
            if(mean(Size_2(:)) == 1)
                Current_Data_Length = 1;
            else
                Current_Data_Length = Size_2(1);
            end
        else
            if(mean(Size_2(:)) == 1)
                Current_Data_Length = Size_1(1);
            else
                error("MCPL_To_MAT : 2D variables found within partition file, expected one dimensional data.");
            end
        end
        Data_Length(Current_Reference) = Current_Data_Length;
    end
    %Check number of partitions from each file matches each other
    if(range(Number_Of_Partitions) ~= 0)
        error("MCPL_To_MAT : Mismatch between Datastore_Partition_Information partitions and the number of partition files.");
    end
    %Check partitions are sequential
    [Partition_Index, Sorted_Index] = sort(Partition_Index, 'ascend');
    if(sum(diff(Partition_Index)==1) ~= numel(Partition_Index) - 1)
        error("MCPL_To_MAT : Sorted Partitions aren't sequential, missing partition files.");
    end
    %Sort Files into order
    Files = Files(Sorted_Index);
    Data_Length = Data_Length(Sorted_Index);
    %% Recreate File_Chunks structure
    File_Chunks = struct('Chunk', num2cell(1:1:length(Files)),'Temp_File_Path', fullfile({Files.folder}, {Files.name}), 'Start', 0, 'End', 0, 'Events', num2cell(Data_Length));
end