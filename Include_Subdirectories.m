%Adds a list of directories used for storing functions in an organised manner to the MATLAB path enviroment
%:Inputs
% - Directory_List (cell array / string / char)
% - Relative_Directory (integer / boolean) ; looks at the calling function(s) relative directories to change root directory. Integer specifies the exact calling function from the tree, boolean uses the parent function directory.
%:Outputs
% - Valid (boolean) ; Whether the function input was valid
function Valid_Input = Include_Subdirectories(Directory_List, Parent_Function_Relative_Directory)
    %% If no input is supplied, add all default subdirectories to path
    switch nargin
        case 0
            %use default list
            disp("No Directories supplied, using default list.");
            Directory_List = {'Data_Operations','File_Operations','Input_Validation','Parpool', 'Waitbar', 'WinRAR'};
            Using_Parent_Function_Relative_Directory = false;
        case 1
            %use user input
            Using_Parent_Function_Relative_Directory = false;
        case 2
            %switch relative directories
            if(isnumeric(Parent_Function_Relative_Directory))
                Using_Parent_Function_Relative_Directory = true;
            elseif(islogical(Parent_Function_Relative_Directory))
                if(Parent_Function_Relative_Directory)
                    Using_Parent_Function_Relative_Directory = true;
                    Parent_Function_Relative_Directory = 2;
                else
                    Using_Parent_Function_Relative_Directory = false;
                end
            end
        otherwise
            %do nothing
            disp("Unhandled Input");
    end
    %Ensure flag for relative directory is set
    if(~exist('Using_Parent_Function_Relative_Directory','var'))
        Using_Parent_Function_Relative_Directory = false;
    end
    %Ensure numeric value is set for relative directories
    if(~exist('Parent_Function_Relative_Directory','var'))
        Parent_Function_Relative_Directory = 0;
    end
    %Verify function input
    if(isstring(Directory_List))
        %convert to cell
        Directory_List = cellstr(Directory_List);
    elseif(ischar(Directory_List))
        %convert to cell
        Directory_List = cellstr(Directory_List);
    else
        %do nothing
    end
    %if the input is now formatted as a cell
    if(iscell(Directory_List))
        Valid_Input = 1;
    else
        %do nothing
        Valid_Input = 0;
    end
    %if input directory list is valid
    if(Valid_Input)
        %Get current root directory
        Current_Root_Directory = pwd;
        if(Using_Parent_Function_Relative_Directory)
            %get calling function tree
            Calling_Functions = dbstack('-completenames');
            if(Parent_Function_Relative_Directory >= length(Calling_Functions))
                Parent_Function_Relative_Directory = length(Calling_Functions);
            end
            [pathstr,~,~] = fileparts(Calling_Functions(Parent_Function_Relative_Directory).file);
        else
            %change directory to matlab script's directory
            [pathstr,~,~] = fileparts(mfilename('fullpath'));
        end
        cd(pathstr);
        clear pathstr name ext;

        %add relative MATLAB Code paths for functions in subdirectory
        for Current_Directory_Number = 1:length(Directory_List)
            %get current directory
            Current_Directory = strcat(cd,filesep,Directory_List{Current_Directory_Number});
            %if directory exists, add to path
            if(isdir(Current_Directory))
                addpath(Current_Directory)
            else
                disp(strcat("Invalid Directory: ",Current_Directory));
            end
        end
        %ensure current directory is also added to path
        addpath(Current_Root_Directory);
        cd(Current_Root_Directory);
    else
        disp("Invalid function input, expected array of string, char or cell datatype");
    end
    %% Output validation
    if nargout == 0
        clear Valid_Input;
    end
end