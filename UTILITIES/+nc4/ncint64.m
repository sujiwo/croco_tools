classdef ncint64 < nc4.nctype
    methods
        function self = ncint64(varargin)
            self@nc4.nctype('NC_INT64', varargin);
        end
    end
end
        
