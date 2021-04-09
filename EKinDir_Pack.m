%% Pack EKinDir vectors from Dx, Dy, Dz and Energy
function [EKinDir_1, EKinDir_2, EKinDir_3] = EKinDir_Pack(Dx, Dy, Dz, Energy)
    %% Input Handling
    if(nargin ~= 4)
        error("EKinDir_Pack : Expected four inputs.");
    end
    if(length(Dx) ~= length(Dy))
        error("EKinDir_Pack : Dx and Dy length mismatch.");
    end
    if(length(Dx) ~= length(Dz))
        error("EKinDir_Pack : Dx and Dz length mismatch.");
    end
    if(length(Dx) ~= length(Energy))
        error("EKinDir_Pack : Dx and Energy length mismatch.");
    end
    %Ensure Dx, Dy, Dz are unit vectors for the packing algorithm (if already unit vectors; this
    %won't change the value of each vector)
    RMS = sqrt(Dx.^2 + Dy.^2 + Dz.^2);
    Dx = Dx./RMS;
    Dy = Dy./RMS;
    Dz = Dz./RMS;
    %output (x,y,sign(z)) - default
    EKinDir_1 = Dx;
    EKinDir_2 = Dy;
    EKinDir_3 = Dz;
    Abs_X = abs(Dx);
    Abs_Y = abs(Dy);
    Condition_1 = abs(Dz) < max(Abs_X, Abs_Y);
    %Invert Z where appropriate conditions are met
    Inv_Z(length(Dz)) = Inf;
    Inv_Z(Condition_1 & ~Floating_Point_Equal(Dz(:),0)) = 1./Dz(Condition_1 & ~Floating_Point_Equal(Dz(:),0));
    Condition_2 = Abs_X >= Abs_Y;
    %output (1/z,y,sign(x))
    EKinDir_1(Condition_1 & Condition_2) = Inv_Z(Condition_1 & Condition_2);
    EKinDir_2(Condition_1 & Condition_2) = Dy(Condition_1 & Condition_2);
    EKinDir_3(Condition_1 & Condition_2) = Dx(Condition_1 & Condition_2);
    %output (x,1/z,sign(y))
    EKinDir_1(Condition_1 & ~Condition_2) = Dx(Condition_1 & ~Condition_2);
    EKinDir_2(Condition_1 & ~Condition_2) = Inv_Z(Condition_1 & ~Condition_2);
    EKinDir_3(Condition_1 & ~Condition_2) = Dy(Condition_1 & ~Condition_2);
    %Encode the sign bit into the energy value
    EKinDir_3(:) = sign(EKinDir_3(:)) .* Energy(:);
end