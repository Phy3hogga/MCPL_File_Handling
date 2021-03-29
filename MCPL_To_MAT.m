%% Reads MCPL file and converts it to a MAT file
function MAT_File_Path = MCPL_To_MAT(MCPL_File_Path, Read_Parameters)
    %% Input handling
    %Argument handling
    if(nargin == 0)
        error("MCPL_To_MAT requires a file path input");
    elseif(nargin == 1)
        warning("Using default settings for MCPL_To_Mat");
    end
    %Skipping decompression
    [Struct_Var_Value, Struct_Var_Valid] = Verify_Structure_Input(Read_Parameters, 'Skip_Uncompress', false);
    if(Struct_Var_Valid)
        Skip_Uncompress = Struct_Var_Value;
    else
        Skip_Uncompress = false;
    end
    %Sorting events by weighting
    [Struct_Var_Value, Struct_Var_Valid] = Verify_Structure_Input(Read_Parameters, 'Sort_Events_By_Weight', true);
    if(Struct_Var_Valid)
        Sort_Events_By_Weight = Struct_Var_Value;
    else
        Sort_Events_By_Weight = true;
    end
    %Removal of 0 weight events
    [Struct_Var_Value, Struct_Var_Valid] = Verify_Structure_Input(Read_Parameters, 'Remove_Zero_Weights', true);
    if(Struct_Var_Valid)
        Remove_Zero_Weights = Struct_Var_Value;
    else
        Remove_Zero_Weights = true;
    end
    %Removal of 0 weight events
    [Struct_Var_Value, Struct_Var_Valid] = Verify_Structure_Input(Read_Parameters, 'Remove_Temp_Files', true);
    if(Struct_Var_Valid)
        Remove_Temp_Files = Struct_Var_Value;
    else
        Remove_Temp_Files = true;
    end
    [Struct_Var_Value, Struct_Var_Valid] = Verify_Structure_Input(Read_Parameters, 'Parpool_Num_Cores', 1);
    if(Struct_Var_Valid)
        Parpool_Num_Cores = Struct_Var_Value;
    else
        Parpool_Num_Cores = 1;
    end
    %% Uncompress GZ archive using WinRAR
    [Directory_Path, Filename, Extension] = fileparts(MCPL_File_Path);
    if(strcmpi(Extension, '.gz'))
        %Extraction of GZ archive (if the file format matches)
        Uncompressed_File_Path = strcat(Directory_Path, filesep, Filename, '-UNCOMPRESS');
        %Only use RAR_Parameters field if it has been parsed by previous settings structure
        if(isfield(Read_Parameters, 'RAR_Parameters'))
            RAR_Parameters = Read_Parameters.RAR_Parameters;
        else
            RAR_Parameters = struct();
        end
        if(~Skip_Uncompress)
            Successful_Uncompress = UNRAR(MCPL_File_Path, Uncompressed_File_Path, RAR_Parameters);
        else
            Successful_Uncompress = 1;
        end
        clear Skip_Uncompress RAR_Parameters;
        if(~Successful_Uncompress)
            error('Error uncompressing GZ !bofile format');
        end
        MCPL_File_List = {};
        File_Path_Search = Uncompressed_File_Path;
        clear Successful_Uncompress Uncompressed_File_Path;
    elseif(strcmpi(Extension, '.mcpl'))
        %No extraction required
        if(isfile(MCPL_File_Path))
            MCPL_File_List{1} = MCPL_File_Path;
        else
            MCPL_File_List = {};
        end
        File_Path_Search = MCPL_File_Path;
    else
        error('Unexpected file format');
    end
    clear Extension Filename Directory_Path MCPL_File_Path;

    %% Find MCPL file(s) that aren't explicitly stated in the original input (if a directory is specified)
    if(isempty(MCPL_File_List))
        File_List = Search_Files(File_Path_Search, '.mcpl');
        if(isempty(fieldnames(File_List)))
            error('No mcpl files found');
        end
    end
    clear File_Path_Search MCPL_File_List;
    
    %Preallocate the list of file paths
    MAT_File_Path{1:length(File_List)} = '';
    %% Read files
    for Read_Index = 1:length(File_List)
        %Path and reference to file
        File_Path = fullfile(File_List(Read_Index).folder, filesep, File_List(Read_Index).name);
        File_ID = fopen(File_Path, 'r');
        %% Read file header
        Header = MCPL_Read_Header(File_ID);
        Endian = Header.Endian;
        %Get size of file
        fseek(File_ID, 0, 'eof');
        File.End = ftell(File_ID);
        File.Data = File.End - Header.End;
        fseek(File_ID, Header.End, 'bof');
        Successful_Close = fclose(File_ID);
        if(Successful_Close == -1)
            warning(["MCPL file failed to close: ", File_List(Read_Index).name]);
        end
        %Find size of a single photon's information from header information
        if(Header.Valid)
            %% Prep for reading data
            %Variable type data
            if(Header.Opt_SinglePrecision)
                Byte_Size = 4;
                Byte_Type = 'single';
            else
                Byte_Size = 8;
                Byte_Type = 'double';
            end
            %Add bit depth information to header information
            Header.Byte_Size = Byte_Size;
            Header.Byte_Type = Byte_Type;
            %Tracking the byte position of specific data within a single photon's data
            Byte_Position = 0;
            Dynamic_Photon_Field = 0;
            %Tracking how many bytes to read for a specific photon in total
            Photon_Byte_Count = 1;
            if(Header.Opt_Polarisation)
                Dynamic_Photon_Field = Dynamic_Photon_Field + 3;
                Photon_Byte_Count = Photon_Byte_Count + 3 * Byte_Size;
                %Polarisation data within the photon byte string
                Byte_Position = Byte_Position + 1;
                Byte_Split.Px.Start = Byte_Position;
                Byte_Position = Byte_Position + Byte_Size - 1;
                Byte_Split.Px.End = Byte_Position;

                Byte_Position = Byte_Position + 1;
                Byte_Split.Py.Start = Byte_Position;
                Byte_Position = Byte_Position + Byte_Size - 1;
                Byte_Split.Py.End = Byte_Position;

                Byte_Position = Byte_Position + 1;
                Byte_Split.Pz.Start = Byte_Position;
                Byte_Position = Byte_Position + Byte_Size - 1;
                Byte_Split.Pz.End = Byte_Position;
            end
            %Addition of extra field; ekindir decompression
            Dynamic_Photon_Field = Dynamic_Photon_Field + 8;
            Photon_Byte_Count = Photon_Byte_Count + 7 * Byte_Size;
            %Position data within the photon byte string
            Byte_Position = Byte_Position + 1;
            Byte_Split.X.Start = Byte_Position;
            Byte_Position = Byte_Position + Byte_Size - 1;
            Byte_Split.X.End = Byte_Position;

            Byte_Position = Byte_Position + 1;
            Byte_Split.Y.Start = Byte_Position;
            Byte_Position = Byte_Position + Byte_Size - 1;
            Byte_Split.Y.End = Byte_Position;

            Byte_Position = Byte_Position + 1;
            Byte_Split.Z.Start = Byte_Position;
            Byte_Position = Byte_Position + Byte_Size - 1;
            Byte_Split.Z.End = Byte_Position;

            %Displacement/Energy Vectors within the photon byte string
            Byte_Position = Byte_Position + 1;
            Byte_Split.EKinDir_1.Start = Byte_Position;
            Byte_Position = Byte_Position + Byte_Size - 1;
            Byte_Split.EKinDir_1.End = Byte_Position;

            Byte_Position = Byte_Position + 1;
            Byte_Split.EKinDir_2.Start = Byte_Position;
            Byte_Position = Byte_Position + Byte_Size - 1;
            Byte_Split.EKinDir_2.End = Byte_Position;

            Byte_Position = Byte_Position + 1;
            Byte_Split.EKinDir_3.Start = Byte_Position;
            Byte_Position = Byte_Position + Byte_Size - 1;
            Byte_Split.EKinDir_3.End = Byte_Position;

            %Time within the photon byte string
            Byte_Position = Byte_Position + 1;
            Byte_Split.Time.Start = Byte_Position;
            Byte_Position = Byte_Position + Byte_Size - 1;
            Byte_Split.Time.End = Byte_Position;

            if(~Header.Opt_UniversalWeight)
                Dynamic_Photon_Field = Dynamic_Photon_Field + 1;
                Photon_Byte_Count = Photon_Byte_Count + Byte_Size - 1;
                %Weight within the photon byte string
                Byte_Position = Byte_Position + 1;
                Byte_Split.Weight.Start = Byte_Position;
                Byte_Position = Byte_Position + Byte_Size - 1;
                Byte_Split.Weight.End = Byte_Position;
            end
            %Fixed size data
            %If PDGCode varies
            if(Header.Opt_UniversalPDGCode == 0)
                Photon_Byte_Count = Photon_Byte_Count + 4;
                Byte_Position = Byte_Position + 1;
                Byte_Split.PDGCode.Start = Byte_Position;
                Byte_Position = Byte_Position + 4 - 1;
                Byte_Split.PDGCode.End = Byte_Position;
            end
            %If userflags are present
            if(Header.Opt_Userflag)
                Photon_Byte_Count = Photon_Byte_Count + 4;
                Byte_Position = Byte_Position + 1;
                Byte_Split.UserFlag.Start = Byte_Position;
                Byte_Position = Byte_Position + 4 - 1;
                Byte_Split.UserFlag.End = Byte_Position;
            end
            %Add additional fields to header
            Header.Photon_Byte_Count = Photon_Byte_Count;
            Header.Byte_Split = Byte_Split;
            Header.Sort_Events_By_Weight = Sort_Events_By_Weight;
            Header.Remove_Zero_Weights = Remove_Zero_Weights;

            %% Create root output directory
            Temp_Output_File_Root = fullfile(File_List(Read_Index).folder, filesep, 'TEMPORARY');
            Directory_Creation_Success = Attempt_Directory_Creation(Temp_Output_File_Root);
            if(~Directory_Creation_Success)
                warning(strcat("Failed to create temporary output directory for: ", File_List(Read_Index).name));
            end
            %% Parallel core processing setup
            if(Parpool_Num_Cores > 1)
                Parpool = Parpool_Create(Parpool_Num_Cores);
            end
            %If parpool is disabled; requires the number of cores to be assigned to a variable
            if(~exist('Parpool', 'var'))
                Parpool.NumWorkers = Parpool_Num_Cores;
            end
            [~, System_Memory] = memory;
            %Calculate file chunk interval depending on available system memory and the total file size
            %Addition of 3 fields for unpacking of EKinDir(E,x,y,z) to E, Dir(x,y,z)
            %Memory compensation factor of 40% use for additional data handling overhead while processing
            Interval_Memory = floor(((System_Memory.PhysicalMemory.Available * 0.4) / Parpool.NumWorkers) / (Photon_Byte_Count + (3 * Byte_Size)));
            Interval_File = floor((File.Data / Parpool.NumWorkers) / (Photon_Byte_Count + (3 * Byte_Size)));
            %Change chunk interval based on memory available and file size
            if(Interval_Memory < Interval_File)
                Interval = Interval_Memory;
            else
                Interval = Interval_File;
            end
            %% Photon Data
            Chunks = 1:Interval:Header.Particles;
            %Edit final chunk (should be minor) to add any remaining photon chunks that aren't included via equal division
            %Either adds an additional chunk or appends a few extra events to the final chunk depending on discrepency
            if(Chunks(end) ~= Header.Particles)
                if(Header.Particles - Chunks(end) > Parpool_Num_Cores)
                    Chunks(end + 1) = Header.Particles;
                else
                    Chunks(end) = Header.Particles;
                end
            end
            %Calculate dynamic and corrected interval
            Interval = Chunks(2:end) - Chunks(1:end-1);
            File_Chunks = struct('Chunk', num2cell(1:1:length(Chunks)-1),'Temp_File_Path', fullfile(strcat(Temp_Output_File_Root, filesep, arrayfun(@num2str, 1:1:length(Chunks)-1, 'UniformOutput', 0), '.mat')), 'Start', num2cell(((Chunks(1:end-1)-1) * Photon_Byte_Count) + Header.End), 'End', num2cell(((Chunks(1:end-1)-1) + Interval - 1) * Photon_Byte_Count + Header.End + 1), 'Events', num2cell(Interval));
            %End of file correction (should be a single Event)
            if(File_Chunks(end).End ~= File.End)
                %Adjust final chunk end if required
                File_Chunks(end).End = File.End;
                %Adjust chunk size as per end of file
                File_Chunks(end).Events = (File_Chunks(end).End - File_Chunks(end).Start)/Photon_Byte_Count;
            end
            %Add chunks to header
            Header.File_Chunks = File_Chunks;
            %Save header
            Header_File_Path = fullfile(Temp_Output_File_Root, 'Header.mat');
            save(Header_File_Path, '-v7.3', '-struct', 'Header');
            %% Read file chunks aand dump them to disk, sorted individual chunks by weighting
            disp("Reading MCPL file");
            if(Parpool_Num_Cores > 1)
                %Parallel processing
                parfor Current_File_Chunk = 1:length(File_Chunks)
                    MCPL_Dump_Data_Chunk(Header, File_Path, File_Chunks(Current_File_Chunk));
                end
                Parpool_Delete();
            else
                %Single core processing
                for Current_File_Chunk = 1:length(File_Chunks)
                    MCPL_Dump_Data_Chunk(Header, File_Path, File_Chunks(Current_File_Chunk));
                end
            end
            
            %% Combine all output file(s)
            disp("Processing MCPL file into MAT file");
            MAT_File_Path{Read_Index} = MCPL_Merge_Chunks(Header, File_Path);

            %% Cleanup temporary files (including all files within)
            if(Remove_Temp_Files)
                Temporary_Files_Removed = rmdir(Temp_Output_File_Root, 's');
            end
        else
            error(strcat("MCPL file format not found for file: ", File_List(Read_Index).name));
        end
    end
end

%% Read MCPL Header
function Header = MCPL_Read_Header(File_ID)
    %Read file header
    File_Header = fread(File_ID, 8, '*char');
    %Verify valid MCPL file
    if(strcmpi(File_Header(1), 'M') && strcmpi(File_Header(2), 'C') && strcmpi(File_Header(3), 'P') && strcmpi(File_Header(4), 'L'))
        %Verify file version compatibility
        Format_Version = (File_Header(5)-'0')*100 + (File_Header(6)-'0')*10 + (File_Header(7)-'0');
        if(~any(Format_Version == [2,3]))
            warning(["MCPL file version may not be compatible with this read script: ", File_List(File_Index).name]);
        end
        %Warning for version 2 files
        if(Format_Version == 2)
            warning(["MCPL file version 2 not fully tested for this script, legacy input: ", File_List(File_Index).name]);
        end
        %Get computer native endian type
        [~, ~, Computer_Endian] = computer;
        %Check file endian type
        File_Endian_Version = File_Header(8);
        Header.Endian_Switch = false;
        %Compare endianness between the file and computer
        if(strcmpi(File_Endian_Version, 'B'))
            if(strcmpi(Computer_Endian, 'L'))
                Header.Endian_Switch = true;
            end
            Endian.T32 = 'b';
            Endian.T64 = 's';
        elseif(strcmpi(File_Endian_Version, 'L'))
            if(strcmpi(Computer_Endian, 'B'))
                Header.Endian_Switch = true;
            end
            Endian.T32 = 'l';
            Endian.T64 = 'a';
        else
            error(["Could not determine endianness for file: ", File_List(File_Index).name]);
        end
        Header.Endian = Endian;
        % Cleanup
        clear File_Header Endian_Version;

        %% Start reading remainder of file content
        %% Numeric Content
        %Number of particles
        Header.MCPL_Version = Format_Version;
        clear Format_Version;
        Header.Particles = fread(File_ID, 1, 'uint64', 0, Endian.T64);
        %Number of Comments
        N_Comments = fread(File_ID, 1, 'uint32', 0, Endian.T32);
        %Number of Blobs
        N_Blobs = fread(File_ID, 1, 'uint32', 0, Endian.T32);
        %Optional Header Information
        Header.Opt_Userflag = fread(File_ID, 1, 'uint32', 0, Endian.T32);
        Header.Opt_Polarisation = fread(File_ID, 1, 'uint32', 0, Endian.T32);
        Header.Opt_SinglePrecision = fread(File_ID, 1, 'uint32', 0, Endian.T32);
        Header.Opt_UniversalPDGCode = fread(File_ID, 1, 'int32', 0, Endian.T32);
        Header.Opt_ParticleSize = fread(File_ID, 1, 'uint32', 0, Endian.T32);
        Header.Opt_UniversalWeight = fread(File_ID, 1, 'uint32', 0, Endian.T32);
        if(Header.Opt_UniversalWeight)
            Header.Opt_UniversalWeightValue = fread(File_ID, 1, 'uint64', 0, Endian.T64);
        end
        %Unspecified Use
        Header.Opt_Signature = 0 + 1*Header.Opt_SinglePrecision + 2*Header.Opt_Polarisation + 4*Header.Opt_UniversalPDGCode + 8*Header.Opt_UniversalWeight + 16*Header.Opt_Userflag;

        %% String Content
        % Data Source
        Header.Source{1} = MCPL_Read_String(File_ID, Endian);
        % Comments
        if(N_Comments > 0)
            for Current_Comment = 1:N_Comments
                Header.Comments{Current_Comment} = MCPL_Read_String(File_ID, Endian);
            end
        else
            Header.Comments{1} = '';
        end

        %Blobs
        if(N_Blobs > 0)
            %Blob Keys
            for Current_Blob = 1:N_Blobs
                Header.Blobs.Key{Current_Blob} = MCPL_Read_String(File_ID, Endian);
            end
            %Blob Data
            for Current_Blob = 1:N_Blobs
                Header.Blobs.Data{Current_Blob} = MCPL_Read_String(File_ID, Endian);
                %Transpose the data stored to a readable format
                Header.Blobs.Data{Current_Blob} = Header.Blobs.Data{Current_Blob}';
            end
        else
            Header.Blobs.Key{1} = '';
            Header.Blobs.Data{1} = '';
        end
        Header.Valid = true;
    else
        Header.Valid = false;
    end
    Header.End = ftell(File_ID);
end

%% Read an MCPL String
function MCPL_String = MCPL_Read_String(File_ID, Endian)
    %Get string length
    String_Length = fread(File_ID, 1, 'uint32', 0, Endian.T32);
    %Read string if length is positive
    if(String_Length > 0)
        MCPL_String = fread(File_ID, String_Length, '*char');
    else
        MCPL_String = '';
    end
end

%% Read and Dump an MCPL File Chunk to MAT file
function MCPL_Dump_Data_Chunk(Header, File_Path, File_Chunk)
    %% Preallocate arrays
    if(Header.Opt_Polarisation)
        Px = zeros(File_Chunk.Events, 1, Header.Byte_Type);
        Py = zeros(File_Chunk.Events, 1, Header.Byte_Type);
        Pz = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    end
    X = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    Y = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    Z = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    Dx = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    Dy = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    Dz = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    Energy = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    Time = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    EKinDir_1 = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    EKinDir_2 = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    EKinDir_3 = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    if(~Header.Opt_UniversalWeight)
        Weight = zeros(File_Chunk.Events, 1, Header.Byte_Type);
    end
    if(Header.Opt_UniversalPDGCode == 0)
        PDGCode = zeros(File_Chunk.Events, 1, 'int32');
    end
    if(Header.Opt_Userflag)
        UserFlag = zeros(File_Chunk.Events, 1, 'uint32');
    end

    %% Read chunk of data from file
    File_ID = fopen(File_Path, 'r');
    %fseek(File_ID, Header.End, 'bof');
    fseek(File_ID, File_Chunk.Start, 'bof');
    Byte_String = fread(File_ID, Header.Photon_Byte_Count * File_Chunk.Events, 'uint8=>uint8');
    Successful_Chunk_Close = fclose(File_ID);
    if(Successful_Chunk_Close == -1)
        warning(strcat("Unsuccessful File Close when Reading Data for Chunk: ", num2str(File_Chunk.Chunk)));
    end
    %% Input Handling
    for Event_Number = 1:File_Chunk.Events
        %Offset for read start
        Photon_Offset = (Event_Number - 1) * Header.Photon_Byte_Count;

        %Variable Data
        if(Header.Opt_Polarisation)
            Px(Event_Number) = typecast(Byte_String(Header.Byte_Split.Px.Start + Photon_Offset : Header.Byte_Split.Px.End + Photon_Offset), Header.Byte_Type);
            Py(Event_Number) = typecast(Byte_String(Header.Byte_Split.Py.Start + Photon_Offset : Header.Byte_Split.Py.End + Photon_Offset), Header.Byte_Type);
            Pz(Event_Number) = typecast(Byte_String(Header.Byte_Split.Pz.Start + Photon_Offset : Header.Byte_Split.Pz.End + Photon_Offset), Header.Byte_Type);
        end
        X(Event_Number) = typecast(Byte_String(Header.Byte_Split.X.Start + Photon_Offset : Header.Byte_Split.X.End + Photon_Offset), Header.Byte_Type);
        Y(Event_Number) = typecast(Byte_String(Header.Byte_Split.Y.Start + Photon_Offset : Header.Byte_Split.Y.End + Photon_Offset), Header.Byte_Type);
        Z(Event_Number) = typecast(Byte_String(Header.Byte_Split.Z.Start + Photon_Offset : Header.Byte_Split.Z.End + Photon_Offset), Header.Byte_Type);
        EKinDir_1(Event_Number) = typecast(Byte_String(Header.Byte_Split.EKinDir_1.Start + Photon_Offset : Header.Byte_Split.EKinDir_1.End + Photon_Offset), Header.Byte_Type);
        EKinDir_2(Event_Number) = typecast(Byte_String(Header.Byte_Split.EKinDir_2.Start + Photon_Offset : Header.Byte_Split.EKinDir_2.End + Photon_Offset), Header.Byte_Type);
        EKinDir_3(Event_Number) = typecast(Byte_String(Header.Byte_Split.EKinDir_3.Start + Photon_Offset : Header.Byte_Split.EKinDir_3.End + Photon_Offset), Header.Byte_Type);
        Time(Event_Number) = typecast(Byte_String(Header.Byte_Split.Time.Start + Photon_Offset : Header.Byte_Split.Time.End + Photon_Offset), Header.Byte_Type);
        if(~Header.Opt_UniversalWeight)
            Weight(Event_Number) = typecast(Byte_String(Header.Byte_Split.Weight.Start + Photon_Offset : Header.Byte_Split.Weight.End + Photon_Offset), Header.Byte_Type);
        end
        %Fixed size data
        if(Header.Opt_UniversalPDGCode == 0)
            PDGCode(Event_Number) = typecast(Byte_String(Header.Byte_Split.PDGCode.Start + Photon_Offset : Header.Byte_Split.PDGCode.End + Photon_Offset), 'int32');
        end
        if(Header.Opt_Userflag)
            UserFlag(Event_Number) = typecast(Byte_String(Header.Byte_Split.UserFlag.Start + Photon_Offset : Header.Byte_Split.UserFlag.End + Photon_Offset), 'int32');
        end
        %% TODO: Adjust endian-ness of byte order
        if(Header.Endian_Switch)
            disp("WARNING: Verify results, system endianness changed");
        end
    end
    %% Unpacking of energy and direction vectors (method depends on MCPL File Version)
    if(Header.MCPL_Version == 3)
        %Unpack energy
        Energy = abs(EKinDir_3);
        sign_bit = sign(EKinDir_3);
        %input is (x,y,sign(z))
        Dx = EKinDir_1;
        Dy = EKinDir_2;
        Dz = sign_bit .* sqrt(max(0.0, 1.0 - (EKinDir_1.^2 + Dy.^2)));

        %input is (x,1/z,sign(y))
        Y_Vec = find(abs(EKinDir_2) > 1.0);
        if(~isempty(Y_Vec))
            Dx(Y_Vec) = EKinDir_1(Y_Vec);
            Dz(Y_Vec) = 1.0 ./ EKinDir_2(Y_Vec);
            Dy(Y_Vec) = sign_bit(Y_Vec) .* sqrt(max(0.0, 1.0 - (EKinDir_1(Y_Vec).^2 + Dz(Y_Vec).^2)));
        end
        %input is (1/z,y,sign(x))
        X_Vec = find(abs(EKinDir_1) > 1.0);
        if(~isempty(X_Vec))
            Dy(X_Vec) = EKinDir_2(X_Vec);
            Dz(X_Vec) = 1.0 ./ EKinDir_1(X_Vec);
            Dx(X_Vec) = sign_bit(X_Vec) .* sqrt(max(0.0,1.0 - (EKinDir_2(X_Vec).^2 + Dz(X_Vec).^2)));
        end
    elseif(Header.MCPL_Version == 2)
        %% UNTESTED; ONLY SOURCE OF DATA FROM LEGACY MCPL COMPONENT VERSION
        Dz = 1.0 - abs(EKinDir_1) - abs(EKinDir_2);
        %Lower Hemisphere
        Hemisphere_Lower = find(Dz < 0);
        if(any(Hemisphere_Lower))
            X_Vec = -ones(size(EKinDir_1),'int8');
            X_Vec(EKinDir_1 >= 0.0) = 1;
            Dx(Hemisphere_Lower) = (1.0 - abs(EKinDir_2(Hemisphere_Lower))) .* X_Vec(Hemisphere_Lower);
            Y_Vec = -ones(size(EKinDir_1),'int8');
            Y_Vec(EKinDir_2 >= 0.0) = 1;
            Dy(Hemisphere_Lower) = (1.0 - abs(EKinDir_1(Hemisphere_Lower))) .* Y_Vec(Hemisphere_Lower);
        end
        %Upper Hemisphere
        Hemisphere_Upper = find(Dz >= 0);
        if(any(Hemisphere_Upper))
            Dx(~Hemisphere) = EKinDir_1(~Hemisphere);
            Dy(~Hemisphere) = EKinDir_2(~Hemisphere);
        end
        %Project from octahedron to unit sphere
        N = 1.0 ./ sqrt(Dx.^2 + Dy.^2 + Dz.^2);
        Dx = N .* Dx;
        Dy = N .* Dy;
        Dz = N .* Dz;
        %Get Energy
        Energy = EKinDir_3;
        %Get sign bit
        Sign_Bit = EKinDir_3(2) < 0;
        Energy(Sign_Bit) = -Energy(Sign_Bit);
        Dy(Sign_Bit) = 0.0;
    else
        warning(strcat("Could not unpack direction and energy for chunk: ", num2str(File_Chunk.Chunk)));
        Dx(:) = NaN;
        Dy(:) = NaN;
        Dz(:) = NaN;
        Energy(:) = NaN;
    end
    %Convert photon energy to KeV
    Energy = Energy ./ 1e-3;

    %% Sort events by weighting
    if(~Header.Opt_UniversalWeight && Header.Sort_Events_By_Weight)
        [Weight, Sorted_Index] = sort(Weight, 'descend');
        Energy = Energy(Sorted_Index);
        Time = Time(Sorted_Index);
        X = X(Sorted_Index);
        Y = Y(Sorted_Index);
        Z = Z(Sorted_Index);
        Dx = Dx(Sorted_Index);
        Dy = Dy(Sorted_Index);
        Dz = Dz(Sorted_Index);
        EKinDir_1 = EKinDir_1(Sorted_Index);
        EKinDir_2 = EKinDir_2(Sorted_Index);
        EKinDir_3 = EKinDir_3(Sorted_Index);
        if(Header.Opt_Polarisation)
            Px = Px(Sorted_Index);
            Py = Px(Sorted_Index);
            Pz = Px(Sorted_Index);
        end
        if(Header.Opt_UniversalPDGCode == 0)
            PDGCode = PDGCode(Sorted_Index);
        end
        if(Header.Opt_Userflag)
            UserFlag = UserFlag(Sorted_Index);
        end
    end

    %% Save file chunk to a temporary file for combination later
    save(File_Chunk.Temp_File_Path, '-v7.3', 'Weight', 'Energy', 'Time', 'X', 'Y', 'Z', 'Dx', 'Dy', 'Dz', 'EKinDir_1', 'EKinDir_2', 'EKinDir_3');
    if(Header.Opt_Polarisation)
        save(File_Chunk.Temp_File_Path, '-v7.3', '-append', 'Px', 'Py', 'Pz');
    end
    if(Header.Opt_UniversalPDGCode == 0)
        save(File_Chunk.Temp_File_Path, '-v7.3', '-append', 'PDGCode');
    end
    if(Header.Opt_Userflag)
        save(File_Chunk.Temp_File_Path, '-v7.3', '-append', 'UserFlag');
    end
end

%% Merge MCPL File Chunks into a single MAT file
function Merged_File_Path = MCPL_Merge_Chunks(Header, File_Path)
    %File_Chunk data loading from header
    File_Chunks = Header.File_Chunks;
    Chunk_Files = {File_Chunks(:).Temp_File_Path};
    %Load matfile references to each file containing a chunk of processed data
    Chunk_Matfile_References = cellfun(@matfile, Chunk_Files, 'UniformOutput', false);
    %Find limits of array sizes for each chunk
    Chunk_Events = [File_Chunks(:).Events];
    Total_Chunk_Events = sum(Chunk_Events(:));
    Max_Chunk_Events = max(Chunk_Events);

    %For each file, split the weights being read into chunks
    Weight_Chunk = floor(Max_Chunk_Events / length(Chunk_Matfile_References));
    
    %Create table for sorting weights
    Table_Fields = {};
    Table_Datatypes = {};
    Table_Fields{end + 1} = 'File_Index';
    Table_Datatypes{end + 1} = 'int64';
    Table_Fields{end + 1} = 'File_Row';
    Table_Datatypes{end + 1} = 'int64';
    if(Header.Opt_Polarisation)
        Table_Fields{end + 1} = 'Px';
        Table_Datatypes{end + 1} = Header.Byte_Type;
        Table_Fields{end + 1} = 'Py';
        Table_Datatypes{end + 1} = Header.Byte_Type;
        Table_Fields{end + 1} = 'Px';
        Table_Datatypes{end + 1} = Header.Byte_Type;
    end
    Table_Fields{end + 1} = 'X';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'Y';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'Z';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'Dx';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'Dy';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'Dz';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'Energy';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'Time';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'EKinDir_1';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'EKinDir_2';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    Table_Fields{end + 1} = 'EKinDir_3';
    Table_Datatypes{end + 1} = Header.Byte_Type;
    if(~Header.Opt_UniversalWeight)
        Table_Fields{end + 1} = 'Weight';
        Table_Datatypes{end + 1} = Header.Byte_Type;
    end
    if(Header.Opt_UniversalPDGCode == 0)
        Table_Fields{end + 1} = 'PDGCode';
        Table_Datatypes{end + 1} = 'int32';
    end
    if(Header.Opt_Userflag)
        Table_Fields{end + 1} = 'UserFlag';
        Table_Datatypes{end + 1} = 'uint32';
    end
    Weight_Table = table('Size', [(Weight_Chunk + 1) * length(Chunk_Matfile_References), length(Table_Fields)],'VariableTypes', Table_Datatypes);
    Weight_Table.Properties.VariableNames = Table_Fields;
    %Pre-fill table with NaN values
    Weight_Table.File_Index(:) = NaN;
    Weight_Table.File_Row(:) = NaN;
    if(Header.Opt_Polarisation)
        Weight_Table.Px(:) = NaN;
        Weight_Table.Py(:) = NaN;
        Weight_Table.Pz(:) = NaN;
    end
    Weight_Table.X(:) = NaN;
    Weight_Table.Y(:) = NaN;
    Weight_Table.Z(:) = NaN;
    Weight_Table.Dx(:) = NaN;
    Weight_Table.Dy(:) = NaN;
    Weight_Table.Dz(:) = NaN;
    Weight_Table.Energy(:) = NaN;
    Weight_Table.Time(:) = NaN;
    Weight_Table.EKinDir_1(:) = NaN;
    Weight_Table.EKinDir_2(:) = NaN;
    Weight_Table.EKinDir_3(:) = NaN;
    if(~Header.Opt_UniversalWeight)
        Weight_Table.Weight(:) = NaN;
    end
    if(Header.Opt_UniversalPDGCode == 0)
        Weight_Table.PDGCode(:) = NaN;
    end
    if(Header.Opt_Userflag)
        Weight_Table.UserFlag(:) = NaN;
    end
    %Store relative position of each chunk in the weight table
    Weight_Table_Position = 1;
    %Track current row within the file
    Current_File_Row = zeros(size(Chunk_Matfile_References));
    %Track remaining rows after a chunk is sorted and some unsorted data remains
    Current_File_Row_Offset = zeros(size(Chunk_Matfile_References));
    %Track unsorted possitions chunk position after sorting
    Sorted_List_Final_Chunk_Position = NaN(size(Chunk_Matfile_References));
    %Logical index for which file needs reading
    Read_Event = ones(size(Chunk_Events), 'logical');

    %Index of current row within the file to write to
    File_Write_Index = 1;
    Removed_Zero_Count = 0;
    %File path to write the merged data to
    [Root_File_Path, Filename, ~] = fileparts(File_Path);
    Merged_File_Path = fullfile(Root_File_Path, strcat(Filename, ".mat"));
    %Create merged file
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
    Merged_File_Reference.EKinDir_1(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.EKinDir_2(Total_Chunk_Events, 1) = Empty_Byte_Type;
    Merged_File_Reference.EKinDir_3(Total_Chunk_Events, 1) = Empty_Byte_Type;
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
    
    %% Read chunks while some files still have rows left to read
    while(any(Read_Event))
        %Re-check files that still need reading
        Read_Values = find(Read_Event == true);
        %If any files need reading still
        if(any(Read_Values))
            for Read_Chunk_Index = 1:length(Read_Values)
                Current_File = Read_Values(Read_Chunk_Index);
                %Exception for initial reading of not having to add an index
                Start_File_Row = Current_File_Row(Current_File) + 1;
                %Don't need to read this file; all entries processed
                if(Start_File_Row == Chunk_Events(Current_File))
                    Read_Event(Current_File) = 0;
                end
                End_File_Row = min([Chunk_Events(Current_File), Start_File_Row + Weight_Chunk - Current_File_Row_Offset(Current_File)]);
                %Only place new data when valid indicies are presented
                if(End_File_Row >= Start_File_Row)
                    %If finished reading the file; skip reading more data in future (also clip to the number of lines left in the file)
                    if(End_File_Row >= Chunk_Events(Current_File))
                        %Don't need to read this file; all entries processed
                        Read_Event(Current_File) = 0;
                        %Verify end of file row is used; overwrite larger values as sanity check
                        End_File_Row = Chunk_Events(Current_File);
                    end
                    %Add results to table
                    End_Weight_Table_Position = Weight_Table_Position + (End_File_Row - Start_File_Row);
                    Weight_Table.File_Index(Weight_Table_Position:End_Weight_Table_Position) = Current_File;
                    Weight_Table.File_Row(Weight_Table_Position:End_Weight_Table_Position) = Start_File_Row:1:End_File_Row;
                    if(Header.Opt_Polarisation)
                        Weight_Table.Px(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Px(Start_File_Row:End_File_Row, 1);
                        Weight_Table.Py(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Py(Start_File_Row:End_File_Row, 1);
                        Weight_Table.Pz(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Pz(Start_File_Row:End_File_Row, 1);
                    end
                    Weight_Table.X(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.X(Start_File_Row:End_File_Row, 1);
                    Weight_Table.Y(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Y(Start_File_Row:End_File_Row, 1);
                    Weight_Table.Z(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Z(Start_File_Row:End_File_Row, 1);
                    Weight_Table.Dx(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Dx(Start_File_Row:End_File_Row, 1);
                    Weight_Table.Dy(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Dy(Start_File_Row:End_File_Row, 1);
                    Weight_Table.Dz(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Dz(Start_File_Row:End_File_Row, 1);
                    Weight_Table.Energy(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Energy(Start_File_Row:End_File_Row, 1);
                    Weight_Table.Time(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Time(Start_File_Row:End_File_Row, 1);
                    Weight_Table.EKinDir_1(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.EKinDir_1(Start_File_Row:End_File_Row, 1);
                    Weight_Table.EKinDir_2(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.EKinDir_2(Start_File_Row:End_File_Row, 1);
                    Weight_Table.EKinDir_3(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.EKinDir_3(Start_File_Row:End_File_Row, 1);
                    if(~Header.Opt_UniversalWeight)
                        Weight_Table.Weight(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.Weight(Start_File_Row:End_File_Row, 1);
                    end
                    if(Header.Opt_UniversalPDGCode == 0)
                        Weight_Table.PDGCode(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.PDGCode(Start_File_Row:End_File_Row, 1);
                    end
                    if(Header.Opt_Userflag)
                        Weight_Table.UserFlag(Weight_Table_Position:End_Weight_Table_Position) = Chunk_Matfile_References{Current_File}.UserFlag(Start_File_Row:End_File_Row, 1);
                    end
                    %Store current row
                    Current_File_Row(Current_File) = End_File_Row;
                    %Iterate by chunk size for next file
                    Weight_Table_Position = End_Weight_Table_Position + 1;
                end
            end
            %Reset write position in Weight_Table to the first row (circular buffer)
            Weight_Table_Position = 1;
            if(Header.Sort_Events_By_Weight)
                %Remove weight events with exactly 0.0 weighting
                if(Header.Remove_Zero_Weights)
                    Remove_Zero_Indicies = find(Floating_Point_Equal(Weight_Table.Weight, 0.0));
                    if(~isempty(Remove_Zero_Indicies))
                        Removed_Zero_Count = Removed_Zero_Count + length(Remove_Zero_Indicies);
                        if(Header.Opt_Polarisation)
                            Weight_Table.Px(Remove_Zero_Indicies) = NaN;
                            Weight_Table.Py(Remove_Zero_Indicies) = NaN;
                            Weight_Table.Pz(Remove_Zero_Indicies) = NaN;
                        end
                        Weight_Table.X(Remove_Zero_Indicies) = NaN;
                        Weight_Table.Y(Remove_Zero_Indicies) = NaN;
                        Weight_Table.Z(Remove_Zero_Indicies) = NaN;
                        Weight_Table.Dx(Remove_Zero_Indicies) = NaN;
                        Weight_Table.Dy(Remove_Zero_Indicies) = NaN;
                        Weight_Table.Dz(Remove_Zero_Indicies) = NaN;
                        Weight_Table.Energy(Remove_Zero_Indicies) = NaN;
                        Weight_Table.Time(Remove_Zero_Indicies) = NaN;
                        Weight_Table.EKinDir_1(Remove_Zero_Indicies) = NaN;
                        Weight_Table.EKinDir_2(Remove_Zero_Indicies) = NaN;
                        Weight_Table.EKinDir_3(Remove_Zero_Indicies) = NaN;
                        if(~Header.Opt_UniversalWeight)
                            Weight_Table.Weight(Remove_Zero_Indicies) = NaN;
                        end
                        if(Header.Opt_UniversalPDGCode == 0)
                            Weight_Table.PDGCode(Remove_Zero_Indicies) = NaN;
                        end
                        if(Header.Opt_Userflag)
                            Weight_Table.UserFlag(Remove_Zero_Indicies) = NaN;
                        end
                    end
                end
                %Sort by weight followed by a sub sort of file row (if two weights are identical but are out of order) in the file
                [Weight_Table, ~] = sortrows(Weight_Table, {'Weight','File_Row'}, {'descend','descend'}, 'MissingPlacement','first');

                %Find the indicies of NaN elements
                NaN_Elements = find(isnan(Weight_Table.Weight));
                %Find the linear index in the sorted list of the last entry from each file
                Sorted_List_Final_Chunk_Position(:) = NaN;
                for Read_Chunk_Index = 1:length(Read_Values)
                    File_Line_Check = Current_File_Row(Read_Chunk_Index);
                    while isnan(Sorted_List_Final_Chunk_Position(Read_Chunk_Index))
                        %Safety catch to stop accidental overwriting
                        if(isnan(Sorted_List_Final_Chunk_Position(Read_Chunk_Index)))
                            File_Line_Check_Index = find((Weight_Table.File_Index == Read_Chunk_Index) & (Weight_Table.File_Row == File_Line_Check));
                            %Assign value if not NaN
                            if(~isnan(Weight_Table.Weight(File_Line_Check_Index)))
                                Sorted_List_Final_Chunk_Position(Read_Chunk_Index) = File_Line_Check_Index;
                            end
                            File_Line_Check = File_Line_Check - 1;
                        end
                    end
                end
                %Verify not all indicies are NaN
                if(~all(isnan(Sorted_List_Final_Chunk_Position)))
                    %Find final valid entry
                    [Last_Entry_For_File, ~] = nanmin(Sorted_List_Final_Chunk_Position);
                    if(Last_Entry_For_File > 1)
                        %can write from index 1:Last_Entry_For_File - 1; as the final element will need comparing with next chunk
                        Write_Final_Index = Last_Entry_For_File - 1;
                    else
                        Write_Final_Index = Last_Entry_For_File;
                    end

                    %Remove any NaN elements (located in top half of table after sorting) from writing to file
                    if(any(NaN_Elements))
                        Write_Initial_Index = NaN_Elements(end) + 1;
                    else
                        Write_Initial_Index = 1;
                    end

                    %Sense-check linear indicies
                    if(Write_Initial_Index <= Write_Final_Index)
                        %Write to file
                        File_Write_Index_End = File_Write_Index + length(Weight_Table.File_Index(Write_Initial_Index:Write_Final_Index)) - 1;

                        if(Header.Opt_Polarisation)
                            Merged_File_Reference.Px(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.File_Row(Write_Initial_Index:Write_Final_Index);
                            Merged_File_Reference.Py(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.File_Row(Write_Initial_Index:Write_Final_Index);
                            Merged_File_Reference.Pz(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.File_Row(Write_Initial_Index:Write_Final_Index);
                        end
                        Merged_File_Reference.X(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.X(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.Y(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Y(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.Z(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Z(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.Dx(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Dx(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.Dy(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Dy(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.Dz(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Dz(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.Energy(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Energy(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.Time(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Time(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.EKinDir_1(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.EKinDir_1(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.EKinDir_2(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.EKinDir_2(Write_Initial_Index:Write_Final_Index);
                        Merged_File_Reference.EKinDir_3(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.EKinDir_3(Write_Initial_Index:Write_Final_Index);
                        if(~Header.Opt_UniversalWeight)
                            Merged_File_Reference.Weight(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Weight(Write_Initial_Index:Write_Final_Index);
                        end
                        if(Header.Opt_UniversalPDGCode == 0)
                            Merged_File_Reference.PDGCode(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.PDGCode(Write_Initial_Index:Write_Final_Index);
                        end
                        if(Header.Opt_Userflag)
                            Merged_File_Reference.UserFlag(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.UserFlag(Write_Initial_Index:Write_Final_Index);
                        end

                        %Increment for next pass
                        File_Write_Index = File_Write_Index_End + 1;

                        %Reset position for inserting next read chunks (stops duplication of data when reaching the end of files if not all overwritten)
                        Weight_Table.File_Index(1:Write_Final_Index) = NaN;
                        Weight_Table.File_Row(1:Write_Final_Index) = NaN;
                        if(Header.Opt_Polarisation)
                            Weight_Table.Px(1:Write_Final_Index) = NaN;
                            Weight_Table.Py(1:Write_Final_Index) = NaN;
                            Weight_Table.Pz(1:Write_Final_Index) = NaN;
                        end
                        Weight_Table.X(1:Write_Final_Index) = NaN;
                        Weight_Table.Y(1:Write_Final_Index) = NaN;
                        Weight_Table.Z(1:Write_Final_Index) = NaN;
                        Weight_Table.Dx(1:Write_Final_Index) = NaN;
                        Weight_Table.Dy(1:Write_Final_Index) = NaN;
                        Weight_Table.Dz(1:Write_Final_Index) = NaN;
                        Weight_Table.Energy(1:Write_Final_Index) = NaN;
                        Weight_Table.Time(1:Write_Final_Index) = NaN;
                        Weight_Table.EKinDir_1(1:Write_Final_Index) = NaN;
                        Weight_Table.EKinDir_2(1:Write_Final_Index) = NaN;
                        Weight_Table.EKinDir_3(1:Write_Final_Index) = NaN;
                        if(~Header.Opt_UniversalWeight)
                            Weight_Table.Weight(1:Write_Final_Index) = NaN;
                        end
                        if(Header.Opt_UniversalPDGCode == 0)
                            Weight_Table.PDGCode(1:Write_Final_Index) = NaN;
                        end
                        if(Header.Opt_Userflag)
                            Weight_Table.UserFlag(1:Write_Final_Index) = NaN;
                        end
                        
                        %Find all remaining instances of file references to find the read offset for each file
                        for Current_File = 1:length(Current_File_Row_Offset)
                            Current_File_Row_Offset(Current_File) = max(0, length(find(Weight_Table.File_Index == Current_File)));
                        end
                    else
                        disp("Warning: Final Index larger than Initial Index when combining files.");
                    end
                else
                    disp("Warning: Determination of position for the final data from chunk undetermined.");
                end
            else
                %add all results, no sorting needed
            end
        end
    end
    %% Write final contents of table
    %Fully delete the NaN allocated table elements now the table can be destroyed from allocated memory
    NaN_Elements = find(isnan(Weight_Table.X) | isnan(Weight_Table.Y) | isnan(Weight_Table.Z));
    Weight_Table(NaN_Elements,:) = [];
    %Write to file
    File_Write_Index_End = File_Write_Index + height(Weight_Table) - 1;
    if(Header.Opt_Polarisation)
        Merged_File_Reference.Px(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.File_Row(:);
        Merged_File_Reference.Py(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.File_Row(:);
        Merged_File_Reference.Pz(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.File_Row(:);
    end
    Merged_File_Reference.X(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.X(:);
    Merged_File_Reference.Y(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Y(:);
    Merged_File_Reference.Z(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Z(:);
    Merged_File_Reference.Dx(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Dx(:);
    Merged_File_Reference.Dy(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Dy(:);
    Merged_File_Reference.Dz(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Dz(:);
    Merged_File_Reference.Energy(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Energy(:);
    Merged_File_Reference.Time(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Time(:);
    Merged_File_Reference.EKinDir_1(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.EKinDir_1(:);
    Merged_File_Reference.EKinDir_2(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.EKinDir_2(:);
    Merged_File_Reference.EKinDir_3(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.EKinDir_3(:);
    if(~Header.Opt_UniversalWeight)
        Merged_File_Reference.Weight(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.Weight(:);
    end
    if(Header.Opt_UniversalPDGCode == 0)
        Merged_File_Reference.PDGCode(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.PDGCode(:);
    end
    if(Header.Opt_Userflag)
        Merged_File_Reference.UserFlag(File_Write_Index:File_Write_Index_End, 1) = Weight_Table.UserFlag(:);
    end

    % Remove additional 0'd entries if they exist
    if(Header.Particles > File_Write_Index_End)
        if(Header.Opt_Polarisation)
            Merged_File_Reference.Px(File_Write_Index_End:Header.Particles, 1) = [];
            Merged_File_Reference.Py(File_Write_Index_End:Header.Particles, 1) = [];
            Merged_File_Reference.Pz(File_Write_Index_End:Header.Particles, 1) = [];
        end
        Merged_File_Reference.X(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.Y(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.Z(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.Dx(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.Dy(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.Dz(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.Energy(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.Time(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.EKinDir_1(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.EKinDir_2(File_Write_Index_End:Header.Particles) = [];
        Merged_File_Reference.EKinDir_3(File_Write_Index_End:Header.Particles) = [];
        if(~Header.Opt_UniversalWeight)
            Merged_File_Reference.Weight(File_Write_Index_End:Header.Particles) = [];
        end
        if(Header.Opt_UniversalPDGCode == 0)
            Merged_File_Reference.PDGCode(File_Write_Index_End:Header.Particles) = [];
        end
        if(Header.Opt_Userflag)
            Merged_File_Reference.UserFlag(File_Write_Index_End:Header.Particles) = [];
        end
    end
    %Display output file progress
    disp(strcat("Input Events             : ", num2str(Header.Particles)));
    disp(strcat("Retained Events          : ", num2str(File_Write_Index_End)));
    disp(strcat("Removed 0 Weight Events  : ", num2str(Removed_Zero_Count)));
    
    %Edit the number of events in the stored file
    Header.Particles = File_Write_Index_End;
    %% Copy header into the data file last (ensures writing is finished, if misssing file is invalid)
    Merged_File_Reference.Header = Header;
end

%% REFERENCE FOR READING GZIP FROM FILESTREAM; UNUSED BUT POTENITAL UPGRADE IN FUTURE
%% https://www.cs.usfca.edu/~parrt/doc/java/JavaIO-notes.pdf
% File_Str = javaObject('java.io.FileInputStream',File_Path);
% inflatedStr = javaObject('java.util.zip.GZIPInputStream', File_Str);
% charStr     = javaObject('java.io.InputStreamReader', inflatedStr);
% lines       = javaObject('java.io.BufferedReader', charStr);
% currentLine = lines.readLine();