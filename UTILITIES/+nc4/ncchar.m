classdef ncchar < nc4.nctype
    methods
        function self = ncchar(varargin)
            self@nc4.nctype('NC_CHAR', varargin);
        end
    end
end

