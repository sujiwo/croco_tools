classdef ncobject < handle
    properties
    end

    methods
        function self = ncobject(vargin)
        end

        % Subsref dispatcher
        function res = subsref(self, tstruct)
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
                            res = self.(subs)(args{:});
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
                        v = var(self, subs);
                        indices = operator(1).subs;
                        v.set(indices, input);
                        return
                    else
                        error("### Unsupported variable index")
                    end
            end
        end
    end
end
