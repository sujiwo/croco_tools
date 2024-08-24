classdef ncobject < handle
    properties
    end

    methods
        function self = ncobject(vargin)
        end

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
                            res = self.(subs)(tstruct(1).subs{:});
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
                    x = 1;
            end
            x = 1;
        end
    end
end
