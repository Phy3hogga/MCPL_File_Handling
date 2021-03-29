%Turn MAT file back into a MCPL file
function MCPL_File = MAT_To_MCPL(Mat_File_Path)
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
                    if(Number_Of_Events ~= Header.Particles)
                        disp(strcat("Warning: Number of events in header disagrees with the quantity of data in MAT file: ", Mat_File_Path));
                        disp(strcat("Number of Events in Header: ", num2str(Header.Particles)));
                        disp(strcat("Number of Events in Data: ", num2str(Number_Of_Events)));
                        disp(strcat("Number of Events that will be Translated to MCPL: ", num2str(min(Number_Of_Events, Header.Particles))));
                        Number_Of_Events = Header.Particles;
                    end
                    
                    %% Create header content
                    Header_Content = 'MCPL003';
                    %Get computer native endian type
                    [~, ~, Computer_Endian] = computer;
                    %Compare endianness between the file and computer
                    if(strcmpi(Computer_Endian, 'B'))
                        Header_Content(end + 1) = 'B';
                    elseif(strcmpi(Computer_Endian, 'L'))
                        Header_Content(end + 1) = 'L';
                    else
                        error(strcat("Could not determine system endianness to write file : ", Mat_File_Path));
                    end
                    
                    
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
    %placeholder
    MCPL_File = Mat_File_Path;
end