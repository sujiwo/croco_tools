classdef netcdf < handle

  properties
    id = -1;
    filename = '';
    finfo = -1;
    ndims, nvars, dims, vars;
  end       

  methods

    function o = netcdf(pth, mode)
      if (nargin==1)
        mode = 'NOWRITE';
      end
      if (isa(pth, 'nc4.netcdf'))
        o.id = pth;
        o.filename = pth.filename;
      elseif (ischar(pth) || isstring(pth))
        o.id = netcdf.open(pth, mode);
        o.filename = pth;
      end
      % Dimensions & Variables
      o.finfo = ncinfo(pth);
      ndims = size(o.finfo.Dimensions(:));
      nvars = size(o.finfo.Variables(:));
      o.ndims = ndims(1);
      o.nvars = nvars(1);
      for i=1:o.ndims
        o.dims{i} = o.finfo.Dimensions(i).Name;
      end
      for i=1:o.nvars
        o.vars{i} = o.finfo.Variables(i).Name;
      end
    end

    function r = close(self)
        netcdf.close(self.id);
        r=1;
    end

    function ds = dimensions(self)
        ds = self.dims;
    end

    function vars = variables(self)
      vars = self.vars;
    end

    function dm = dim(self, dimName)
      dimId = netcdf.inqDimID(self.id, dimName);
      dm = nc4.ncdim(self.id, dimId);
    end

    function ret = var(self, varName)
      varId = netcdf.inqVarID(self.id, varName);
      ret = nc4.ncvar(self.id, varId);
    end

    function redef(self)
        netcdf.reDef(self.id);
    end

    function endef(self)
        netcdf.endDef(self.id);
    end

    function a = att(self, name)
        a = nc4.ncatt(self.id, netcdf.getConstant("NC_GLOBAL"), name);
    end

    function dm = createDimension(self, name, length)
      dm = nc4.ncdim.create(self.id, name, length);
    end

    function vr = createVariable(self, name, ntype)
      vr = nc4.ncvar.create(self.id, name, ntype);
    end

    function at = createAttribute(self, name, val)
        at = nc4.ncatt.create(self.id, name);
      if nargin==3
        at.set(val);
      end
    end

  end

  methods(Static)
      function fileobj = create(path, permission)
          if nargin==1
              permission=netcdf.getConstant('CLOBBER');
          end
          permission = bitor(permission, netcdf.getConstant('NETCDF4'));
          id = netcdf.create(path, permission);
          netcdf.close(id);
          fileobj = nc4.netcdf(path, 'write');
      end
  end


end
