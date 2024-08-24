classdef ncvar < nc4.ncobject

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
       s@nc4.ncobject();
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
      if l==1, return, end
      order = linspace(l, 1, l);
      ct = permute(ct, order);
    end

%     XXX: Handle multi-dimensional arrays
    function s = set(self, val)
      if (isscalar(val) && prod(self.size())~=1)
        val2 = zeros(self.size(), nc4.nctype.matlab_type(self.xtype));
        val2(:) = val;
        val = val2;
      end
      netcdf.putVar(self.ncid, self.varId, val);
      s = self;
    end

    function sz = size(self)
      sz = zeros(1, length(self.dimIds));
      dms = self.dimensions();
      for i = 1:length(dms)
        sz(i) = dms{i}.len;
      end
      sz = flip(sz);
    end

    function dms = dimensions(self)
        dms = cell(1, length(self.dimIds));
        for i=1:length(self.dimIds)
          dms{i} = nc4.ncdim(self.ncid, self.dimIds(i));
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

