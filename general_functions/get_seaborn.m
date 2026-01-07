function [ palette,palette_str ] = get_seaborn(  )
    range = 1;
    chosen_palette = 'deep';

    switch chosen_palette 
        case 'deep'
            pal = {'#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3','#937860', '#DA8BC3', '#8C8C8C', '#CCB974', '#64B5CD','#000000'};
    end
    palette = cellfun(@(x) hex2rgb(x,range),pal,'Uniformoutput',0);
    palette_str = cellfun(@(x) sprintf('\\color[rgb]{%.4f %.4f %.4f}',hex2rgb(x,range)),pal,'Uniformoutput',0); 


function [ rgb ] = hex2rgb(hex,range)
%% Input checks:
assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 
if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
%% Tweak inputs if necessary: 
if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end
if nargin == 1
    range = 1; 
end
%% Convert from hex to rgb: 
switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
end

end