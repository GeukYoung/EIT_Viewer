function [output] = FxRecon_AT(data,opt)
    [RM, isbnd] = FxRM_AT_v1;
    
    if nargin == 0
        output = size(RM,1);
    elseif nargin == 1
        if ~ischar(data)
            output = RM*data;
        elseif strcmpi(data,'nelem')
            output = size(RM,1);
        elseif strcmpi(data,'isbnd')
            output = isbnd;
        end
    else
        if strcmpi(opt,'nelem')
            output = size(RM,1);
        elseif strcmpi(opt,'image')
            output = RM*data;
        elseif strcmpi(opt,'sum')
            output = -sum(RM)*data;
        end
    end
end

