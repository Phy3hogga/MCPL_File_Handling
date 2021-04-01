%Turn MAT file back into a MCPL file
function MCPL_File_Path = MAT_To_MCPL(Mat_File_Path, MCPL_File_Path)
    %% Input handling
    if(nargin == 1)
        %% Formulate a MCPL file path derived from the MAT file path
        [Directory, Filename, ~] = fileparts(Mat_File_Path);
        Attempted_MCPL_File_Path = strcat(Directory, filesep, strcat(Filename), '.MCPL');
        if(isfile(Attempted_MCPL_File_Path))
            Non_Duplicate_File_Iteration = 1;
            while(isfile(Attempted_MCPL_File_Path))
                Attempted_MCPL_File_Path = strcat(Directory, filesep, strcat(Filename),'-', num2str(Non_Duplicate_File_Iteration), '.MCPL');
                Non_Duplicate_File_Iteration = Non_Duplicate_File_Iteration + 1;
                if(Non_Duplicate_File_Iteration == 1000)
                    error("Could not formulate a unique filename for MCPL file that doesn't exist");
                end
            end
        end
        MCPL_File_Path = Attempted_MCPL_File_Path;
        warning("Output MCPL filepath not specified");
        disp(strcat("Filepath : ", Attempted_MCPL_File_Path));
        clear Directory Filename Attempted_MCPL_File_Path Non_Duplicate_File_Iteration;
    elseif(nargin == 2)
        %Do nothing; both inputs fulfilled
    else
        error("Expected input of MAT file path");
    end
    %Verify the MAT file exists
    if(isfile(Mat_File_Path))
        %Identify variables within the MAT file
        try
            Mat_File_Variables = who('-file', Mat_File_Path);
        catch
            error(strcat("Invalid MAT file format: ", Mat_File_Path));
        end
        %Verify the header exists
        if(ismember({'Header'}, Mat_File_Variables))
            Field_List = {};
            %% Validate File body content
            %Verify the basic MCPL variables exists
            if (all(ismember({'X','Y','Z','Dx', 'Dy', 'Dz','Weight'}, Mat_File_Variables)))
                %Add verified fields to field list
                Field_List = {'X','Y','Z','Dx', 'Dy', 'Dz','Weight'};
                %Verify EKinDir variables exists
                if(all(ismember({'EKinDir_1','EKinDir_2','EKinDir_3'}, Mat_File_Variables)))
                    Field_List(length(Field_List) + 1:length(Field_List) + 3) = {'EKinDir_1','EKinDir_2','EKinDir_3'};
                    Recalculate_EKinDir = false;
                else
                    %Require energy to attempt recalculation of EKinDir components
                    if(ismember({'Energy'}, Mat_File_Variables))
                        Field_List(length(Field_List) + 1) = {'Energy'};
                        Recalculate_EKinDir = true;
                    else
                        error(strcat("Insufficient Data to calculate EKinDir from Energy and vector components for MCPL file: ", Mat_File_Path));
                    end
                end
                %Verify Time variable exists (will fill with 0's if it doesn't exist)
                if(ismember({'Time'}, Mat_File_Variables))
                    Field_List(length(Field_List) + 1) = {'Time'};
                    Replace_Time = false;
                else
                    Replace_Time = true;
                end
                
                %% Verify length of data for all fields matches the header
                Mat_File_Reference = matfile(Mat_File_Path);
                %Get length of each field
                Size_1 = zeros(length(Field_List),1);
                Size_2 = zeros(length(Field_List),1);
                for Current_Field = 1:length(Field_List)
                    [Size_1(Current_Field), Size_2(Current_Field)] = size(Mat_File_Reference, Field_List{Current_Field});
                end
                %If the number of events in each variable is static
                if (std(Size_1(:)) == 0 && std(Size_2(:)) == 0)
                    %Get the number of events used
                    Number_Of_Events = max(min(Size_1(:)),min(Size_2(:)));
                    clear Size_1 Size_2;
                    
                    %% Load the file header
                    Header = load(Mat_File_Path, '-mat', 'Header');
                    Header = Header.Header;
                    %Verify the number of events matches the length of the data
                    if(Number_Of_Events ~= Header.Particles)
                        disp(strcat("Warning: Number of events in header disagrees with the quantity of data in MAT file: ", Mat_File_Path));
                        disp(strcat("Number of Events in Header: ", num2str(Header.Particles)));
                        disp(strcat("Number of Events in Data: ", num2str(Number_Of_Events)));
                        disp(strcat("Number of Events that will be Translated to MCPL: ", num2str(min(Number_Of_Events, Header.Particles))));
                        Number_Of_Events = Header.Particles;
                    end
                    
                    %% Create header content
                    %% Open file for writing
                    File_ID = fopen(MCPL_File_Path, 'w');
                    %Write file header
                    disp("Writing MCPL Header");
                    File_ID = MCPL_Write_Header(File_ID, Header);
                    %Find memory limits for translating file contents in memory between datatypes
                    [~, System_Memory] = memory;
                    Interval = floor((System_Memory.PhysicalMemory.Available * 0.35) / (Header.Photon_Byte_Count + (3 * Header.Byte_Size)));
                    Chunks = 1:Interval:Header.Particles;
                    if(length(Chunks) > 1)
                        %Edit final chunk (should be minor) to add any remaining photon chunks that aren't included via equal division
                        %Either adds an additional chunk or appends a few extra events to the final chunk depending on discrepency
                        if(Chunks(end) ~= Header.Particles)
                            Chunks(end) = Header.Particles;
                        end
                        %Calculate dynamic and corrected interval
                        Interval = Chunks(2:end) - Chunks(1:end-1);
                        File_Chunks = struct('Chunk', num2cell(1:1:length(Chunks)-1), 'Start', num2cell(((Chunks(1:end-1)-1) * Header.Photon_Byte_Count) + Header.End), 'End', num2cell(((Chunks(1:end-1)-1) + Interval - 1) * Header.Photon_Byte_Count + Header.End + 1), 'Events', num2cell(Interval));
                        %End of file correction (should be a single Event)
                        if(File_Chunks(end).End ~= File.End)
                            %Adjust final chunk end if required
                            File_Chunks(end).End = File.End;
                            %Adjust chunk size as per end of file
                            File_Chunks(end).Events = (File_Chunks(end).End - File_Chunks(end).Start)/Photon_Byte_Count;
                        end
                    else
                        File_Chunks(1).Chunk = 1;
                        File_Chunks(1).Start = 1;
                        File_Chunks(1).End = Header.Particles;
                        File_Chunks(1).Events = Header.Particles;
                    end
                    Max_Chunk_Events = max([File_Chunks(:).Events]);
                    %Preallocate variables
                    if(Header.Opt_SinglePrecision)
                        Empty_Byte_Type = single(0);
                    else
                        Empty_Byte_Type = double(0);
                    end
                    if(Header.Opt_Polarisation)
                        Px(Max_Chunk_Events, 1) = Empty_Byte_Type;
                        Py(Max_Chunk_Events, 1) = Empty_Byte_Type;
                        Pz(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    end
                    X(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    Y(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    Z(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    Dx(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    Dy(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    Dz(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    Energy(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    Time(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    EKinDir_1(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    EKinDir_2(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    EKinDir_3(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    if(~Header.Opt_UniversalWeight)
                        Weight(Max_Chunk_Events, 1) = Empty_Byte_Type;
                    end
                    if(Header.Opt_UniversalPDGCode == 0)
                        PDGCode(Max_Chunk_Events, 1) = int32(0);
                    end
                    if(Header.Opt_Userflag)
                        UserFlag(Max_Chunk_Events, 1) = uint32(0);
                    end
                    %% Read / Write from the MAT to MCPL file in chunks
                    %Reading from MAT file
                    disp("Writing MCPL Data");
                    for Current_File_Chunk = 1:length(File_Chunks)
                        %Read file chunk
                        disp(strcat("Loading MAT File Chunk ", num2str(Current_File_Chunk), " / ", num2str(length(File_Chunks))));
                        if(Header.Opt_Polarisation)
                            Px(:) = Mat_File_Reference.Px(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                            Py(:) = Mat_File_Reference.Py(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                            Pz(:) = Mat_File_Reference.Pz(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        end
                        %Convert m to cm
                        X(:) = Mat_File_Reference.X(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1).*100;
                        Y(:) = Mat_File_Reference.Y(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1).*100;
                        Z(:) = Mat_File_Reference.Z(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1).*100;
                        Dx(:) = Mat_File_Reference.Dx(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        Dy(:) = Mat_File_Reference.Dy(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        Dz(:) = Mat_File_Reference.Dz(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        %If replacing time with 0 values
                        if(Replace_Time)
                            disp("Replacing Time with 0 Values");
                            Time(:) = 0;
                        else
                            Time(:) = Mat_File_Reference.Time(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        end
                        Dz(:) = Mat_File_Reference.Dz(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        %Only need to directly read energy if recalculating EKinDir
                        if(Recalculate_EKinDir)
                            disp("Recalculating EKinDir from Energy and Direction Vectors");
                            %Read Energy in KeV, translate into MeV
                            Energy(:) = Mat_File_Reference.Energy(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1) .* 1e-3;
                            %Unpack EKinDir
                            [EKinDir_1, EKinDir_2, EKinDir_3] = EKinDir_Pack(Dx, Dy, Dz, Energy);
                        else
                            EKinDir_1(:) = Mat_File_Reference.EKinDir_1(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                            EKinDir_2(:) = Mat_File_Reference.EKinDir_2(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                            EKinDir_3(:) = Mat_File_Reference.EKinDir_3(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        end
                        if(~Header.Opt_UniversalWeight)
                            Weight(:) = Mat_File_Reference.Weight(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        end
                        if(Header.Opt_UniversalPDGCode == 0)
                            PDGCode(:) = Mat_File_Reference.PDGCode(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        end
                        if(Header.Opt_Userflag)
                            UserFlag(:) = Mat_File_Reference.UserFlag(File_Chunks(Current_File_Chunk).Start:File_Chunks(Current_File_Chunk).End, 1);
                        end
                        %% Write Data
                        Progress_Steps = 10;
                        Progress_Value = 100*linspace(0, 1, Progress_Steps + 1);
                        Progress_Count = 1;
                        %Display starting progress
                        disp(strcat("Chunk Write Progress : ", num2str(Progress_Value(1) * Progress_Steps), "%"));
                        for Current_Line = 1:length(X)
                            if(Header.Opt_Polarisation)
                                fwrite(File_ID, Px(Current_Line), Header.Byte_Type);
                                fwrite(File_ID, Py(Current_Line), Header.Byte_Type);
                                fwrite(File_ID, Pz(Current_Line), Header.Byte_Type);
                            end
                            fwrite(File_ID, X(Current_Line), Header.Byte_Type);
                            fwrite(File_ID, Y(Current_Line), Header.Byte_Type);
                            fwrite(File_ID, Z(Current_Line), Header.Byte_Type);
                            fwrite(File_ID, EKinDir_1(Current_Line), Header.Byte_Type);
                            fwrite(File_ID, EKinDir_2(Current_Line), Header.Byte_Type);
                            fwrite(File_ID, EKinDir_3(Current_Line), Header.Byte_Type);
                            fwrite(File_ID, Time(Current_Line), Header.Byte_Type);
                            if(~Header.Opt_UniversalWeight)
                                fwrite(File_ID, Weight(Current_Line), Header.Byte_Type);
                            end
                            if(Header.Opt_UniversalPDGCode == 0)
                                fwrite(File_ID, PDGCode(Current_Line), 'int32');
                            end
                            if(Header.Opt_Userflag)
                                fwrite(File_ID, UserFlag(Current_Line), 'uint32');
                            end
                            %Display progress
                            if(mod(Current_Line, round(length(X)/Progress_Steps)) == 0)
                                Progress_Count = Progress_Count + 1;
                                disp(strcat("Chunk Write Progress : ", num2str(Progress_Value(round(Current_Line/length(X)*Progress_Steps) + 1)), "%"));
                            end
                        end
                        %Catch case that final index doesn't get called (due to rounding errors)
                        if(Progress_Count ~= length(Progress_Value))
                            disp(strcat("Chunk Write Progress : ", num2str(Progress_Value(end), "%"));
                        end
                    end
                    %% Close file for writing
                    fclose(File_ID);
                else
                    error(strcat("Unclear data correspondance between variables (different length) unable to compile MCPL file: ", Mat_File_Path));
                end
            else
                error(strcat("Missing Data from MAT file to compile MCPL file: ", Mat_File_Path));
            end
        else
            error(strcat("Missing Header from MAT file to compile MCPL file: ", Mat_File_Path));
        end
    else
        error(strcat("Could not find the specified MAT file: ", Mat_File_Path));
    end
end

%% Write MCPL Header
function File_ID = MCPL_Write_Header(File_ID, Header)
    %Get computer native endian type
    [~, ~, Computer_Endian] = computer;
    %Compare endianness between the file and computer
    if(strcmpi(Computer_Endian, 'B'))
        Endian_Char = 'B';
    elseif(strcmpi(Computer_Endian, 'L'))
        Endian_Char = 'L';
    else
        error(strcat("Could not determine system endianness to write file : ", Mat_File_Path));
    end
    %Number of comments and blobs
    N_Comments = length(Header.Comments);
    N_Blobs = length(Header.Blobs.Key);
    %Write the main header content
    fwrite(File_ID, strcat('MCPL003', Endian_Char), 'char*1');
    fwrite(File_ID, Header.Particles, 'uint64');
    fwrite(File_ID, N_Comments, 'uint32');
    fwrite(File_ID, N_Blobs, 'uint32');
    fwrite(File_ID, Header.Opt_Userflag, 'uint32');
    fwrite(File_ID, Header.Opt_Polarisation, 'uint32');
    fwrite(File_ID, Header.Opt_SinglePrecision, 'uint32');
    fwrite(File_ID, Header.Opt_UniversalPDGCode, 'uint32');
    fwrite(File_ID, Header.Opt_ParticleSize, 'uint32');
    fwrite(File_ID, Header.Opt_UniversalWeight, 'uint32');
    if(Header.Opt_UniversalWeight)
        fwrite(File_ID, Header.Opt_UniversalWeightValue, 'uint64');
    end
    File_ID = MCPL_Write_String(File_ID, Header.Source{1});
    % Comments
    if(N_Comments > 0)
        for Current_Comment = 1:N_Comments
            File_ID = MCPL_Write_String(File_ID, Header.Comments{Current_Comment});
        end
    end
    %Blobs
    if(N_Blobs > 0)
        %Blob Keys
        for Current_Blob = 1:N_Blobs
            File_ID = MCPL_Write_String(File_ID, Header.Blobs.Key{Current_Blob});
        end
        %Blob Data
        for Current_Blob = 1:N_Blobs
            File_ID = MCPL_Write_String(File_ID, Header.Blobs.Data{Current_Blob});
        end
    end
end

%% Write MCPL string
function File_ID = MCPL_Write_String(File_ID, String)
    String_Length = length(String);
    fwrite(File_ID, String_Length, 'uint32');
    fwrite(File_ID, String, '*char');
end