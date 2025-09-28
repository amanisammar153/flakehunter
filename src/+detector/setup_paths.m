function setup_paths()
%SETUP_PATHS Adds project folders to MATLAB path using relative locations.
    root = fh_root();
    addpath(root);
    addpath(genpath(fullfile(root, "+detector")));
    addpath(genpath(fullfile(root, "+tools")));
    addpath(genpath(fullfile(root, "Parameters")));
    addpath(genpath(fullfile(root, "tests")));
    addpath(genpath(fullfile(root, "utils")));
end
