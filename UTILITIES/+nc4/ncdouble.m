classdef ncdouble < nc4.nctype
    methods
        function self = ncdouble(varargin)
            self@nc4.nctype('NC_DOUBLE', varargin);
        end
    end
end

                
