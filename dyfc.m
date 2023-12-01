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
            TC                     % a vector of the fMRI time-coouse
            cW                     % conidtion weights 
            names                  % ROI names
            wTC                    % Windowed TC
            dfc                    % dynamic functional connectivity
            beta                   % beta values first level analysis 
            WinSize (1,1) double   % an integer determining the size of sliding window 
            Overlap (1,1) double   % a percentage of overlaping windows
    end

    methods
        function obj = dyfc(TC, names, cW, varargin)
            obj.TC = TC;
            obj.cW = cW;
            obj.names = names;
            if size(varargin) == 1
                obj.WinSize = varargin{1};
                 obj.Overlap = .5;
            elseif size(varargin) == 2
                obj.WinSize = varargin{1};
                obj.Overlap = varargin{2};
            elseif size(varargin)>2 
                error('Too many input varaibles')
            else
                obj.WinSize = 5;
                obj.Overlap = .5;
            end
        end % function constructor 

        function obj = connectivity(obj) 
            tmp = obj.TC;
            cWtpm = obj.cW;
           if  mod(size(tmp, 1), obj.WinSize)  ~= 0
                fprintf('padding data with zereos to accoumduate the win size.\n')
                nxi = (floor(size(tmp, 1)/obj.WinSize)+1)*obj.WinSize;
                npads = nxi - size(tmp, 1);
                tmp = [zeros(floor(npads/2),size(tmp, 2)); tmp ; zeros(ceil(npads/2),size(tmp, 2))]; 
                cWtpm = [zeros(floor(npads/2),size(obj.cW, 2)); obj.cW ; zeros(ceil(npads/2),size(obj.cW, 2))]; 
           end % if mod
            
           tmp1 = reshape(tmp, [obj.WinSize, size(tmp, 1)/obj.WinSize, size(tmp,2)]);
           tmp2 = [tmp(1+floor(obj.WinSize*obj.Overlap):end,:); zeros(floor(obj.WinSize*obj.Overlap),size(tmp,2))];
           tmp2 = reshape(tmp2, [obj.WinSize, size(tmp, 1)/obj.WinSize, size(tmp2,2)]);


           cWtpm1 = reshape(cWtpm, [obj.WinSize, size(cWtpm, 1)/obj.WinSize, size(cWtpm,2)]);
           cWtpm2 = [cWtpm(1+floor(obj.WinSize*obj.Overlap):end,:); zeros(floor(obj.WinSize*obj.Overlap),size(cWtpm,2))];
           cWtpm2 = reshape(cWtpm2, [obj.WinSize, size(cWtpm, 1)/obj.WinSize, size(cWtpm2,2)]);


           obj.wTC = repmat(hann(obj.WinSize), 1, 2*size(tmp, 1)/obj.WinSize, size(tmp,2)).*reshape(cat(1, tmp1, tmp2), obj.WinSize, 2*size(tmp, 1)/obj.WinSize, size(tmp,2));
           
           obj.cW = reshape(cat(1, cWtpm1, cWtpm2 ), obj.WinSize, 2*size(cWtpm, 1)/obj.WinSize, size(cWtpm,2));
          
           clear dFC

           w = max(0,mean(repmat(hann(obj.WinSize), 1, size(obj.cW, 2)).*obj.cW));
           idx = find( w>0 );

            for itime = 1:size(obj.wTC,2)
                   obj.dfc(:,:, itime) = corr(squeeze(obj.wTC(:,itime,:)));
           end % for itme
                

           data_reduced = obj.dfc(:,:,idx);
           w_reduced = w(idx);
           
           for i1 = 1:size(data_reduced,1)
               for i2 = i1:size(data_reduced,2)
                    obj.beta(i1,i2) = [w_reduced]'\squeeze(data_reduced(i1,i2,:));
               end
           end
             obj.beta =  .5*(obj.beta + obj.beta');     
                      % condition specific
        end
    end

end 