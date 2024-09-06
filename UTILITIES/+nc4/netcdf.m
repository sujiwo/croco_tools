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
                    
                    elseif ismethod(self, subs)
                        % Matlab allows a class method to be called like a.method
                        % (without parenthesis).
                        % We don't allow that
                        args = input(1).subs;
                        if length(args)>0
                            res = self.(subs)(args{:});
                        else
                            if strcmp(subs,'redef')||strcmp(subs,'endef')||strcmp(subs,'close')
                                self.(subs)();
                                return;
                            else
                                res = self.(subs)();
                            end
                        end
                        if (length(input)>1 && isa(res,'nc4.ncobject'))
                            % Shift arguments
                            input(1) = [];
                            varargout{:} = subsref(res, input);
                        else
                            varargout{:} = res;
                        end
                    else
                        varargout{:} = att(self, subs).get();
                        return
                    end
                case '{}'
                    subs = subs{1};
                    v = var(self, subs);
                    if isempty(input)
                        varargout{:} = v; return, end
                    varargout{:} = subsref(v, input);
                case '()'
                    subs = subs{1};
                    d = dim(self, subs);
                    if isempty(input)
                        varargout{:} = d; return, end
                    varargout = subsref(d, input);
                otherwise
                    error("### Unsupported indexing")
            end
        end

        function res = subsasgn(self, operator, input)
            s = operator(1);
            type = s.type;
            subs = s.subs;
            operator(1) = [];

            switch type
                case '.'
                    if isprop(self, subs)
                        self.(subs) = input;
                    elseif ismethod(self, subs)
                        args = operator(1).subs;
                        res = self.(subs)(args{:});
                        if (isa(res, 'nc4.ncobject'))
                            operator(1)=[];
                            subsasgn(res, operator, input);
                        end
                    else
                        try
                            res = att(self, subs).get();
                        catch me
                            % we let on-the-fly attribute creation
                            res = nc4.ncatt.create(self.id, subs);
                        end
                        res.set(input);
                    end
                case '{}'
                    v = self.var(subs{1});
                    subsasgn(v, operator, input);
                case '()'
                    dm = self.dim(subs{1});
                    subsasgn(dm, operator, input);
                otherwise
                    error("Unsupported assignment")
            end

            res = self;
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
