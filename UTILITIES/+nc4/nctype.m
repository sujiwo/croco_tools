classdef nctype

    properties
        type
        dimIds
        dimensions
    end

    methods
        function self = nctype(type, varargin)
            self.type = type;
            if isempty(varargin)
                error("dimension input is empty");
            end
            self.dimensions = varargin{1};
        end
    end

end