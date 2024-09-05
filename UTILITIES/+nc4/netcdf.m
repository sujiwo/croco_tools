classdef netcdf < nc4.ncobject

    properties
        id = -1;
        filename = '';
        finfo = -1;
        ndims, nvars, dims, vars;
    end

    methods

        function o = netcdf(pth, mode)
            o@nc4.ncobject();
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

        function delete(self)
            close(self);
        end

        function close(self)
            netcdf.close(self.id);
            self.id = -1;
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

        function varargout = subsref(self, input)
            s = input(1);
            type = s.type;
            subs = s.subs;
            input(1)=[];

            switch type
                case '.'
                    if isprop(self, subs)
                        varargout = {};
                        varargout{1} = self.(subs);
                        return
                    end
                    if ismethod(self, subs)
                        args = other(1).subs;
                        if length(args)>0
                            varargout{:} = self.(subs)(args{:});
                        else
                            if subs=='redef'|subs=='endef'|subs=='close'
                                self.(subs)();
                            else
                                varargout = self.(subs)();
                            end
                        end
                        return
                    else
                        varargout{:} = att(self, subs).get();
                        return
                    end
                case '{}'
                    subs = subs{1};
                    v = var(self, subs);
                    if isempty(other)
                        varargout{:} = v; return, end
                    varargout = subsref(v, other);
                case '()'
                    subs = subs{1}
                    d = dim(self, subs);
                    if isempty(other)
                        varargout{:} = d; return, end
                    varargout = subsref(d, other);
                otherwise
                    error("### Unsupported indexing")
            end
        end
    end

    methods(Static)
        function fileobj = create(path, permission)
            if nargin==1
                permission=netcdf.getConstant('CLOBBER');
            end
            % Warning: race condition
            delete(path);
            id = netcdf.create(path, permission);
            netcdf.close(id);
            fileobj = nc4.netcdf(path, 'write');
        end
    end


end
