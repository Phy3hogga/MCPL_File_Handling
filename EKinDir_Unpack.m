%% Unpacking of energy and direction vectors (method depends on MCPL File Version)
function [Dx, Dy, Dz, Energy] = EKinDir_Unpack(EKinDir_1, EKinDir_2, EKinDir_3, MCPL_File_Version)
    %% Input Handling
    if(nargin ~=4)
        error("EKinDir_Unpack : Expected four inputs.");
    end
    if(length(EKinDir_1) ~= length(EKinDir_2))
        error("EKinDir_Unpack : EKinDir_1 and EKinDir_2 length mismatch.");
    end
    if(length(EKinDir_1) ~= length(EKinDir_3))
        error("EKinDir_Unpack : EKinDir_1 and EKinDir_3 length mismatch.");
    end
    
    %% Preallocate output with NaN
    Dx(length(EKinDir_1)) = NaN;
    Dy(length(EKinDir_1)) = NaN;
    Dz(length(EKinDir_1)) = NaN;
    Energy(length(EKinDir_1)) = NaN;
    
    %% Unpacking algorithm
    if(MCPL_File_Version == 3)
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
    elseif(MCPL_File_Version == 2)
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
        Sign_Bit = EKinDir_3(:) < 0;
        Energy(Sign_Bit) = -Energy(Sign_Bit);
        Dy(Sign_Bit) = 0.0;
    else
        warning("EKinDir_Unpack : Unknown MCPL File Version.");
    end
end