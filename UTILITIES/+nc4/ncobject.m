classdef ncobject < handle
    properties
    end

    methods
        function self = ncobject(vargin)
        end


        % Subsref dispatcher
        function varargout = subsrefX(self, tstruct)
            s = tstruct(1);
            type = s.type;
            subs = s.subs;
            tstruct(1)=[];

            switch type
                case '.'
                    if ischar(subs)
                        if isprop(self, subs)
                            varargout = {};
                            varargout{1} = self.(subs);
                            return
                        end
                        if ismethod(self, subs)
                            args = tstruct(1).subs;
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
                            varargout = att(self, subs).get();
                            return
                        end
                    end
                case '{}'
                    if ismethod(self, 'var')
                        subs = subs{1};
                        v = var(self, subs);
                        if isempty(tstruct)
                            varargout = v; return, end
                        varargout = subsref(v, tstruct);
                    else
                        error("### Unsupported indexing method")
                    end
                case '()'
                    if ismethod(self, 'get')
                        varargout = get(self, subs{:});
                        return
                    elseif ismethod(self, 'dim')
                        a = dim(self, subs{1});
                        subsref(a, tstruct);
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
