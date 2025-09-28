function r = fh_root()
%FH_ROOT Returns the absolute path of the FlakeHunter project root.
% Works from any subfolder and on any OS.
    persistent cached
    if ~isempty(cached) && isfolder(cached)
        r = cached; return;
    end
    % Strategy:
    % 1) If this file exists, use its location.
    % 2) Else use current folder (last resort).
    here = mfilename('fullpath');
    if ~isempty(here)
        r = fileparts(here);  % this file is in the project root
    else
        r = pwd;
    end
    cached = r;
end
