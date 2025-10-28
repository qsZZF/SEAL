% eegplugin_Newtopoplot - topograph plugin
% Usage:
%   >> eegplugin_csp(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
function vers = eegplugin_Newtopoplot( fig, try_strings, catch_strings)

vers = 'Newtopoplot1.1';
% input check
if nargin < 3
    error('eegplugin_iclabel requires 3 arguments');
end


% create menu
menu = findobj(fig, 'tag', 'plot');

NewtopoMenu = uimenu(menu,'label','New topoplot');
% add new submenu
uimenu( NewtopoMenu, 'label', 'Channel location by name', 'callback', 'pop_Newtopoplot(EEG,1);');
uimenu( NewtopoMenu, 'label', 'Channel location by number', 'callback', 'pop_Newtopoplot(EEG,2);');
uimenu( NewtopoMenu, 'label', 'Topo GUI', 'callback', 'pop_Newtopoplot(EEG);');


end