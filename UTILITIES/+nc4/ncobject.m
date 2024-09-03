classdef ncobject < handle
    properties
    end

    methods
        function self = ncobject(vargin)
        end


        function varargout = subsref(self, other)
            s = other(1);
            type = s.type;
            subs = s.subs;
            other(1) = [];

            switch type
                case '.'
                    if ischar(subs)
                        if isprop(self, subs)
                            res = self.(subs); return, end
                        if ismethod(self, subs)
                            args = tstruct(1).subs;
                            if length(args)>0
                                res = self.(subs)(args{:});
                            else
                                if subs=='redef'|subs=='endef'|subs=='close'
                                    self.(subs)();
                                else
                                    res = self.(subs)();
                                end
                            end
                            return
                        else
                            res = att(self, subs).get();
                            subsref(res, other);
                        end
                    else
                        error("Invalid subscript");
                    end
                case '{}'

                case '()'

                otherwise
                    error ("Unsupported type");
            end
        end

        % Subsref dispatcher
        function res = subsrefX(self, tstruct)
            s = tstruct(1);
            type = s.type;
            subs = s.subs;
            tstruct(1)=[];

            switch type
                case '.'
                    if ischar(subs)
                        if isprop(self, subs)
                            res = self.(subs); return, end
                        if ismethod(self, subs)
                            args = tstruct(1).subs;
                            if length(args)>0
                                res = self.(subs)(args{:});
                            else 
                                if subs=='redef'|subs=='endef'|subs=='close'
                                    self.(subs)();
                                else
                                    res = self.(subs)(); 
                                end
                            end
                            return
                        else
                            res = att(self, subs).get();
                            return
                        end
                    end
                case '{}'
                    if ismethod(self, 'var')
                        subs = subs{1};
                        v = var(self, subs);
                        if isempty(tstruct)
                            res = v; return, end
                        res = subsref(v, tstruct);
                    else
                        error("### Unsupported indexing method")
                    end
                case '()'
                    if ismethod(self, 'get')
                        res = get(self, subs{:});
                        return
                    else
                        error("### Unsupported reference")
                    end
            end
        end

        % subsasgn dispatcher
        function res = subsasgn(self, operator, input)
            s = operator(1);
            type = s.type;
            subs = s.subs;
            operator(1) = [];

            switch type
                case '.'
                    if isprop(self, subs)
                        self.(subs) = input;
                    else
                        try
                            res = att(self, subs).get();
                        catch me
                            % we let on-the-fly attribute creation
                            if isprop(self, 'varId')
                                res = nc4.ncatt.create(self.ncid, subs, self.varId);
                            else
                                res = nc4.ncatt.create(self.id, subs);
                            end
                        end
                        res.set(input);
                    end
                case '{}'
                    if ismethod(self, 'var')
                        subs = subs{1};
                        try
                            v = var(self, subs);
                            % Second-level subsasgn. No recursive
                            switch operator.type
                                case '()'
                                    indices = operator(1).subs;
                                    v.set(indices, input);
                                case '.'
                                    attr = att(v, operator(1).subs);
                                    attr.set(input);
                                otherwise
                            end
                        catch me
                            % Create variable
                            v = self.createVariable(subs, input);
                        end
                    else
                        error("### Unsupported variable index")
                    end
                case '()'
                    if ismethod(self, 'createDimension')
                        att = self.createDimension(subs{1}, input);
                    else
                        error("### Unable to create dimension: not a netCDF object")
                    end
            end
            res = self;
        end
    end
end
