classdef ncvar

  properties
    ncid
    varId
    attributes
    xtype
    dimIds
    numAttrs
    name
    dimNames = {}
  end

  methods
    function s = ncvar(ncId, varId)
       s.ncid = ncId;
       s.varId = varId;
       [s.name, s.xtype, s.dimIds, s.numAttrs] = netcdf.inqVar(ncId, varId);

       for i=0:s.numAttrs-1
         s.attributes{i+1} = netcdf.inqAttName(s.ncid, s.varId, i);
       end
       for i=1:length(s.dimIds)
         [s.dimNames{i}, dlen] = netcdf.inqDim(s.ncid, s.dimIds(i));
       end
    end

    function ct = val(self)
      ct = self.get();
    end

    function ct = get(self)
      ct = netcdf.getVar(self.ncid, self.varId);
      l = length(self.dimIds);
      order = linspace(l, 1, l);
      ct = permute(ct, order);
    end

    function s = set(self, val)
        s = 1;
    end

    function dms = dimensions(self)
        dms = []
        for i=1:length(self.dimIds)
          dms(i) = nc4.ncdim(self.ncid, self.dimIds(i))
        end
    end

    function a = att(self, attrName)
       a = nc4.ncatt(self.ncid, self.varId, attrName);
    end

    function a = createAttribute(self, attrName, val)
% Create and (optionally) set an attribute
      a = nc4.ncatt.create(self.ncid, attrName, self.varId);
      if nargin==3
        a.set(val);
      end
    end

  end

  methods (Static)
% Usage example:
% nc4.ncvar.create(26, 'name', ncdouble('time', 'lat', 'lon'));
%
      function nvar = create(ncid, name, ntype)
% doesn't work in Octave
%          arguments
%              ncid int32
%              name string
%              ntype nc4.nctype
%          end

          ndim = length(ntype.dimensions);
          dimids = zeros(1,ndim);
          for i=1:ndim
              dimids(i) = netcdf.inqDimID(ncid, ntype.dimensions{i});
          end
          varid = netcdf.defVar(ncid, name, ntype.type, dimids);
          nvar = nc4.ncvar(ncid, varid);
      end
  end
end

