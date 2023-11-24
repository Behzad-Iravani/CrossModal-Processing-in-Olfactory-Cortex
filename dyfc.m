% -*- UTF-8 -*-
% The dyfc class conatains the data and method for peroforming sliding
% window dynamic functional analysis
% Copyright (C) Behzad Iravani
% behzadiravani@gmail.com
% 
% Department of Neurology and Neurological sciences, Stanford University, Palo Alto 
% 
% November, 2023 -- Philadelphia
% -------------------------------------------------------------------------
classdef dyfc

    properties
            TC          % a vector of the fMRI time-coouse
            WinSize     % an integer determining the size of sliding window 
            Overlap 
    end

end 