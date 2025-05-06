function h = boxplot2(varargin)
%BOXPLOT2 Enhanced boxplot plots
% 
% h = boxplot2(y)
% h = boxplot2(y,x)
% h = boxplot2(..., p1, v1, ...)
%
% I don't like the original boxplot function... it messes with axis
% dimensions, replaces tick labels with text (even if no rotation is
% needed), and in general hijacks the axis too much.  Plus it doesn't
% return the handles to the plotted objects in an easy-to-use way (I can
% never remember which row corresponds to which part of the boxplot).  This
% version creates more light-handed boxplots, assuming that any cosmetic
% changes (color, tick labels, line specs, etc) can be added by the user
% afterwards if necessary.  It also allows one to create clustered
% boxplots, similar to an unstacked bar graph.
%
% Input variables:
%
%   y:              either a ndata x nx array (as in boxplot) or nx x ny x
%                   ndata array where nx indicates the number of
%                   x-clusters, ny the number of boxes per cluster, and
%                   ndata the number of points per boxplot.
%
%   x:              vector of x locations for each box cluster
%
% Optional input variables (passed as parameter/value pairs)
%
%   notch:          'on' or 'off' ['off']
%
%   orientation:    'vertical' or 'horizontal' ['vertical']
%
%   barwidth:       Barwidth value used to position boxes (see bar) [0.8]
%
%   whisker:        whisker length factor (see boxplot) [1.5]
%
%   axes:           axis to plot to [gca]
%
% Output variables:
%
%   h:              structure of handles to boxplot, with the following
%                   fields: 
%                   'box':      box
%                   'ladj':     lower adjacent value
%                   'lwhis':    lower whisker
%                   'med':      median
%                   'out':      outliers
%                   'uadj':     upper adjacent value
%                   'uwhis':    upper whisker    
%

% Copyright 2012 Kelly Kearney       

% Parse input

p = inputParser;
p.addRequired('y', @isnumeric);
p.addOptional('x', [], @isnumeric);
p.addParamValue('notch', 'off', @ischar);
p.addParamValue('orientation', 'vertical', @ischar);
p.addParamValue('axes', gca, @(x) isscalar(x) && ishandle(x) && strcmp(get(x,'type'),'axes'));
p.addParamValue('barwidth', 0.8, @(x) isscalar(x) && x > 0 && x <= 1);
p.addParamValue('whisker', 1.5, @(x) isscalar(x));

p.parse(varargin{:});

In = p.Results;
In.notch = validatestring(In.notch, {'on', 'off'});
In.orientation = validatestring(In.orientation, {'vertical', 'horizontal'});

if ndims(In.y) == 2
    In.y = permute(In.y, [2 3 1]);
end
[nx, ny, ndata] = size(In.y);

if isempty(In.x)
    In.x = 1:nx;
end

ybox = reshape(In.y, [], ndata)';

% Use bar graph to get x positions

figtmp = figure('visible', 'off');
hax = axes;
hb = bar(In.x, In.y(:,:,1), In.barwidth);

% Version-specific code for getting x positions and boxwidth
if verLessThan('matlab', '8.4.0')
    for ib = 1:length(hb)
        xbar = get(get(hb(ib), 'children'), 'xdata');
        xb(:,ib) = mean(xbar,1);
    end
    boxwidth = diff(minmax(xbar(:,1)));
elseif verLessThan('matlab', '24.1.0')
    for ib = 1:length(hb)
        xb(:,ib) = hb(ib).XData + hb(ib).XOffset;
    end
    boxwidth = diff([hb(1:2).XOffset])*In.barwidth;
else
    % MATLAB 2024a and later
    xb = zeros(nx, ny);
    for ib = 1:length(hb)
        xb(:,ib) = hb(ib).XEndPoints;
    end
    
    if ny > 1
        boxwidth = min(diff(xb(1,:)));
    else
        boxwidth = In.barwidth * 0.8; % Fallback if there's only one box per x position
    end
end



% Create boxplot
hbox = boxplot(ybox, 'positions', xb(:), ...
              'notch', In.notch, ...
              'orientation', In.orientation, ...
              'symbol', '+', ...
              'widths', boxwidth, ...
              'whisker', In.whisker);

% Check if boxplot was created successfully
if isempty(hbox)
    error('boxplot2:NoBoxplotCreated', 'Boxplot creation failed. Check your input data.');
end

% Find and copy objects with detailed error checking
tags = {'Box', 'Lower Adjacent Value', 'Lower Whisker', 'Median', 'Outliers', 'Upper Adjacent Value', 'Upper Whisker'};
fields = {'box', 'ladj', 'lwhis', 'med', 'out', 'uadj', 'uwhis'};

for i = 1:length(tags)
    objects = findall(hax, 'tag', tags{i});
    if isempty(objects)
        warning('boxplot2:NoObjectsFound', 'No objects found with tag "%s".', tags{i});
        h.(fields{i}) = NaN(ny, nx);
    else
        h.(fields{i}) = copyobj(objects, In.axes);
        if isempty(h.(fields{i}))
            warning('boxplot2:CopyFailed', 'Failed to copy objects with tag "%s".', tags{i});
            h.(fields{i}) = NaN(ny, nx);
        else
            h.(fields{i}) = reshape_handles(h.(fields{i}), ny, nx);
        end
    end
end

close(figtmp);

end

function reshaped = reshape_handles(handles, ny, nx)
    num_handles = numel(handles);
    if num_handles == ny * nx
        reshaped = reshape(flipud(handles), ny, nx);
    else
        reshaped = NaN(ny, nx);
        for i = 1:min(num_handles, ny*nx)
            [row, col] = ind2sub([ny, nx], i);
            reshaped(row, col) = handles(i);
        end
        warning('boxplot2:ReshapeFailed', 'Unable to reshape %d handles to %dx%d. Created a padded array.', num_handles, ny, nx);
    end
end

