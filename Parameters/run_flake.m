function out = run_flake(imageDir, varargin)
%RUN_FLAKE  One-click runner for your dataset (no algorithm change).
% Usage:
%   out = run_flake('C:\...\matlab-detector\flakesSET1')
%   out = run_flake('flakesSET1','ConfThresh',0.7,'SizeThreshUm2',600,'RadiusScale',1.15,'OutDir',"results_tuned")

if nargin==0
    imageDir = fullfile(pwd,'flakesSET1');   % default
end

% ===== tuned defaults (UPDATE THESE THREE after your sweep) =====
CONF_TUNED   = 0.70;     % e.g. from sweep T(1,:).ConfThresh
SIZE_TUNED   = 600;      % e.g. from sweep T(1,:).SizeThreshUm2
RADIUS_TUNED = 1.15;     % e.g. from sweep T(1,:).RadiusScale
OUTDIR_DEF   = "results";

% optional overrides
ip = inputParser;
ip.addParameter('ConfThresh',   CONF_TUNED, @(x)isnumeric(x)&&isscalar(x));
ip.addParameter('SizeThreshUm2',SIZE_TUNED, @(x)isnumeric(x)&&isscalar(x));
ip.addParameter('RadiusScale',  RADIUS_TUNED, @(x)isnumeric(x)&&isscalar(x));
ip.addParameter('OutDir',       OUTDIR_DEF, @(s)ischar(s)||isstring(s));
ip.addParameter('SaveOverlays', true, @islogical);
ip.addParameter('SaveMasks',    true, @islogical);
ip.addParameter('SaveCrops',    true, @islogical);
ip.addParameter('Parallel',     false, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

material  = "graphene";   % use the preset that worked best
thickness = "90nm";
MAG       = 20;

% 1) ensure Parameters exist (bootstraps if missing)
if ~isfolder(fullfile(pwd,'Parameters','Flatfields')) || ...
   ~isfolder(fullfile(pwd,'Parameters','GMM_Parameters')) || ...
   ~isfolder(fullfile(pwd,'Parameters','Scan_Magnification'))
    detector.bootstrapParameters(material, thickness, MAG, pwd, imageDir, ...
        'K',3, 'UsedChannels',"BGR", 'MaxSamples',30000, 'RandomSeed',42);
end

% 2) load base parameters (do NOT scale radius here; we'll use RadiusScale knob)
P = detector.loadParams(material, thickness, MAG, ...
        'UsedChannels',"BGR", 'StdThresh',6, ...
        'ConfThresh',opt.ConfThresh, 'SizeThreshUm2',opt.SizeThreshUm2);

% (optional) cache the preset for reproducibility
presetDir = fullfile('Parameters','presets');
if ~exist(presetDir,'dir'), mkdir(presetDir); end
save(fullfile(presetDir,'preset_graphene_90nm_20x.mat'),'P');

% 3) run the batch
out = run_detection_on_dataset(imageDir, ...
    'Params', P, ...
    'OutDir', string(opt.OutDir), ...
    'CleanOutDir', true, ...         % wipe old results first
    'SaveOverlays', opt.SaveOverlays, ...
    'SaveMasks',    opt.SaveMasks, ...
    'SaveCrops',    opt.SaveCrops, ...
    'MaxBackgroundValue', 255, ...
    'RadiusScale',  opt.RadiusScale, ...
    'Parallel',     opt.Parallel, ...
    'Verbose', true);

% Pretty print where things went
fprintf('\nOverlays: %s\n', fullfile(pwd, char(opt.OutDir), 'overlays'));
fprintf('Crops:    %s\n', fullfile(pwd, char(opt.OutDir), 'crops'));
fprintf('Ranking:  %s\n\n', fullfile(pwd, char(opt.OutDir), 'flakes_ranking.csv'));

% 4) auto-evaluate if GT exists
gtCsv = fullfile('gt','gt_boxes.csv');
if isfile(gtCsv)
    R = tools.evaluate_detections(char(opt.OutDir), gtCsv, 0.30);
    try disp(R.perImage); catch, end
end
end
