function [ varargout ] = get_unitinfo( info,all_unitinfo, unitnum )
% function [ varargout ] = get_unitinfo( info,all_unitinfo, unitnum )
%   Get some information about a particular unit or all units.

   
    %disp(all_unitinfo.labels)
    if ~iscell(info)
        info={info};
    end
    N = numel(info);
    if ~exist('unitnum','var') || isempty(unitnum)
       unitnum = NaN; 
    end

    if ~exist('all_unitinfo','var') || isempty(all_unitinfo)
       load all_unitinfo
    end
    
    for i = 1:N
        new_label = info{i};
        colnum = strcmp(new_label,all_unitinfo.labels);
        if sum(colnum)==0
            error('invalid unit label')
        end
        if ~isnan(unitnum)
           varargout(i) = {all_unitinfo.data(unitnum,colnum)};
        else
           varargout(i) = {all_unitinfo.data(:,colnum)};
        end
    end
end
