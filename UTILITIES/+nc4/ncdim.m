classdef ncdim
    properties
        ncid
        dim_id
        name
        len
    end

    methods
        function s = ncdim(ncid, dim_id)
            s.ncid = ncid;
            s.dim_id = dim_id;
            [s.name, s.len] = netcdf.inqDim(ncid, dim_id);
        end

        function sz = size(self)

        end
    end

    methods(Static)
        function dm = create(ncid, name, len)
            if nargin==2 || len==-1
                len = netcdf.getConstant('NC_UNLIMITED');
            end
            dmid = netcdf.defDim(ncid, name, len);
            dm = nc4.ncdim(ncid, dmid);
        end
    end
end
