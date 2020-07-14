% vRF_loadRoot.m
%
% root directory of project files (where *project* folders live)
%
% TCS 7/13/2020

function vRF_root = vRF_loadRoot()


if ~ispc
    thishome = getenv('HOME');
    vRF_root = sprintf('%s/labshare/projects/',thishome);
else
    vRF_root = 'Z:/projects/';
end

return
