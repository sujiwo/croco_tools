classdef ncvar < nc4.ncobject

    properties
        ncid
        varId
        attributes
        xtype
        dimIds
        numAttrs
        name
        myOrientation
        dimNames = {}
        dimSizes
    end

    methods
        function s = ncvar(ncId, varId)
            s@nc4.ncobject();
            s.ncid = ncId;
            s.varId = varId;
            [s.name, s.xtype, s.dimIds, s.numAttrs] = netcdf.inqVar(ncId, varId);
            s.dimIds = fliplr(s.dimIds);

            for i=0:s.numAttrs-1
                s.attributes{i+1} = netcdf.inqAttName(s.ncid, s.varId, i);
            end
            s.dimSizes = zeros(length(s.dimIds),1);
            for i=1:length(s.dimIds)
                [s.dimNames{i}, s.dimSizes(i)] = netcdf.inqDim(s.ncid, s.dimIds(i));
            end
        end

        function ct = val(self)
            ct = self.get();
        end

        function o = orient(self, theOrientation)
            if nargin > 1
                if isempty(theOrientation)
                    theOrientation = 1:length(self.ncsize());
                end
                self.myOrientation = theOrientation;
                o = self;
            else
                o = self.myOrientation;
                if isempty(o)
                    o = 1:length(self.ncsize());
                end
            end
        end

        function ct = getAll(self)
            ct = netcdf.getVar(self.ncid, self.varId);
            l = length(self.dimIds);
            if l==1, return, end
            order = linspace(l, 1, l);
            ct = permute(ct, order);
        end

        function [start,count,stride] = process_indices(self, indices)
            theSize = self.ncsize();
            for i = 1:length(indices)
                if isnumeric(indices{i})
                    if any(diff(diff(indices{i})))
                        disp(' ## Indexing strides must be positive and constant.')
                        return
                    end
                end
            end

            % Flip and permute indices before proceeding,
            %  since we are using virtual indexing.

            theOrientation = self.orient();
            if any(theOrientation < 0) || any(diff(theOrientation) ~= 1)
                for i = 1:length(theOrientation)
                    if theOrientation(i) < 0
                        if isa(indices{i}, 'double')   % Slide the indices.
                            indices{i} = fliplr(theSize(i) + 1 - indices{i});
                        end
                    end
                end
                indices(abs(theOrientation)) = indices;
                theSize(abs(theOrientation)) = theSize;
            end

            if prod(theSize) > 0
                if isempty(theSize)
                    start = 0;
                else
                    start = zeros(1, length(theSize));
                end

                count = ones(1, length(theSize));
                stride = ones(1, length(theSize));
                for i = 1:min(length(indices), length(theSize))
                    k = indices{i};
                    if ~isstr(k) && ~strcmp(k, ':') && ~strcmp(k, '-')
                        start(i) = k(1)-1;
                        count(i) =  length(k);
                        d = 0;
                        if length(k) > 1, d = diff(k); end
                        stride(i) = max(d(1), 1);
                    else
                        count(i) = -1;
                        if i == length(indices) && i < length(theSize)
                            j = i+1:length(theSize);
                            count(j) = -ones(1, length(j));
                        end
                    end
                end
                start(start < 0) = 0;
                stride(stride < 0) = 1;
                for i = 1:length(count)
                    if count(i) == -1
                        maxcount = fix((theSize(i)-start(i)+stride(i)-1) ./ stride(i));
                        count(i) = maxcount;
                    end
                end
                count(count < 0) = 0;
            else
                start = [];
                count = [];
                stride = [];
            end
        end

        % XXX: Should refactor to put start/count/stride out of this function
        function result = get(self, varargin)
            result = [];
            indices = varargin;
            [start,count,stride] = self.process_indices(indices);

            if any(count == 0), error(' ## Bad count.'), end
            if all(count==1)
                start = fliplr(start);
                result = netcdf.getVar(self.ncid, self.varId, start);
                status = 1;
            elseif (all(stride==1))
                start = fliplr(start(:)');
                count = fliplr(count(:)');
                result = netcdf.getVar(self.ncid, self.varId, start, count);
                status = 1;
            else
                start = fliplr(start(:)');
                count = fliplr(count(:)');
                stride = fliplr(stride(:)');
                result = netcdf.getVar(self.ncid, self.varId, start, count, stride);
                status = 1;
            end

            if status>=0 && prod(size(result)) > 0 && (ndims(result)==2) && (strcmp(class(result),'char')) && any(find(size(result)==1))
                %
                % If the read operation was successful
                % and if something was actually returned
                % and if that something has exactly two dimensions
                % and if that something was character
                % and if that character string is actually 1D (ndims never returns 1)
                % then do not permute.
                %
                % This way 1D character arrays are loaded as column vectors.
                %
                % Now if you'll excuse me, after writing this code fragment, I have to go
                % wash my hands vigorously for a few hours (get it off, get it off, get it off, unclean..)
                ;
            elseif status >= 0 && prod(size(result)) > 0
                result = permute(result, length(size(result)):-1:1);
                theOrientation = orient(self);
                if any(theOrientation < 0) | any(diff(theOrientation) ~= 1)
                    for i = 1:length(theOrientation)
                        if theOrientation(i) < 0
                            result = flip(result, abs(theOrientation(i)));
                        end
                    end
                    if length(theOrientation) < 2
                        theOrientation = [theOrientation 2];
                    end
                    result = permute(result, abs(theOrientation));
                end
            elseif status >= 0 && prod(size(result)) == 0
                result = [];
            else
                warning(' ## ncvar/subsref failure.')
            end
        end

        %     XXX: Handle multi-dimensional arrays
        function s = set(self, indices, value)
            [start,count,stride] = self.process_indices(indices);

            if any(count==0), error(" ## Bad count"), end
            while length(count) < 2, count = [count 1]; end
            temp = zeros(count);
            count = count(1:length(start));
            temp(:) = value;
            theOrientation = orient(self);
            if any(theOrientation < 0) || any(diff(theOrientation) ~= 1)
                if length(theOrientation) < 2
                    theOrientation = [theOrientation 2];
                end
                temp = ipermute(temp, abs(theOrientation));
                for i = 1:length(theOrientation)
                    if theOrientation(i) < 0
                        temp = flip(temp, abs(theOrientation(i)));
                    end
                end
            end
            temp = permute(temp, length(size(temp)):-1:1);

            if all(count==1)
                if isempty(start)
                    netcdf.putVar(self.ncid, self.varId, temp);
                else
                    start = fliplr((start(:))');
                    netcdf.putVar(self.ncid, self.varId, start, temp);
                end
            elseif all(stride==1)
                start = fliplr((start(:))');
                count = fliplr((count(:))');
                netcdf.putVar(self.ncid, self.varId, start, count, temp);
            else
                start = fliplr((start(:))');
                count = fliplr((count(:))');
                stride = fliplr((stride(:))');
                netcdf.putVar(self.ncid, self.varId, start, count, stride, temp);
            end
        end

        function sz = ncsize(self)
            sz = self.dimSizes;
        end

        function dt = datatype(self)
            dt = nc4.nctype.matlab_type(self.xtype);
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
                        args = input(1).subs;
                        if ~isempty(args)
                            res = self.(subs)(args{:});
                        else
                            res = self.(subs)();
                        end
                        varargout{:} = res;
                    else
                        varargout{:} = att(self, subs).get();
                        return
                    end
                case '()'
                    varargout{:} = get(self, subs{:});
                otherwise
                    error("Unsupported indexing for variables")
            end
        end
    
        function res = subsasgn(self, operator, input)
            s = operator(1);
            type = s.type;
            subs = s.subs;
            operator(1)=[];
            
            switch type
                % Attribute
                case '.'
                    try
                        a = self.att(subs);
                        a.set(input);
                    catch me
                        self.createAttribute(subs, input);
                    end
                case '()'
                    self.set(subs, input);
                otherwise
                    error("Unsupported mode of assignment")
            end

            res = self;
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

