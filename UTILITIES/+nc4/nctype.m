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

    methods (Static)
        function typ = matlab_type(netcdf_type_id)
            switch netcdf_type_id
                case 1
                    typ = 'int8';
                case 2
                    typ = 'char';
                case 3
                    typ = 'int16';
                case 4
                    typ = 'int32';
                case 5
                    typ = 'single';
                case 6
                    typ = 'double';
                case 7
                    typ = 'uint8';
                case 8
                    typ = 'uint16';
                case 9
                    typ = 'uint32';
                case 10
                    typ = 'int64';
                case 11
                    typ = 'uint64';
                case 12
                    typ = 'string';
            end
        end
    end
end
