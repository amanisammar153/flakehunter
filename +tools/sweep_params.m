function T = sweep_params(imageDir, gtCsv, varargin)
% T = tools.sweep_params('flakesSET1','gt/gt_boxes.csv', ...
%       'Material',"graphene",'Thickness',"90nm",'Magnification',20);
%
% Sweeps small grids over ConfThresh, SizeThreshUm2, RadiusScale.
% Writes each run to results_sweep_<i> and evaluates IoU=0.30.
% Returns a table sorted by F1 desc.

ip = inputParser;
ip.addRequired('imageDir', @(s)ischar(s)||isstring(s));
ip.addRequired('gtCsv',   @(s)ischar(s)||isstring(s));
ip.addParameter('Material',"graphene");
ip.addParameter('Thickness',"90nm");
ip.addParameter('Magnification',20);
ip.addParameter('UsedChannels',"BGR");
ip.addParameter('StdThresh',5);
ip.addParameter('Conf', [0.6 0.7 0.8]);
ip.addParameter('Size', [400 600 800]);
ip.addParameter('Radius',[1.0 1.15 1.3]);
ip.parse(imageDir, gtCsv, varargin{:});
o = ip.Results;

imageDir = char(o.imageDir);
gtCsv    = char(o.gtCsv);

% grid
Conf  = o.Conf(:);
SizeU = o.Size(:);
Rad   = o.Radius(:);

rows = {};
runId = 0;
for ic = 1:numel(Conf)
    for isz = 1:numel(SizeU)
        for ir = 1:numel(Rad)
            runId = runId + 1;
            outDir = sprintf('results_sweep_%02d', runId);
            fprintf('(%d) Conf=%.2f  Size=%d  Radius=%.2f  -> %s\n', ...
                    runId, Conf(ic), SizeU(isz), Rad(ir), outDir);

            % run detection
            out = run_detection_on_dataset(imageDir, ...
                'Material',o.Material,'Thickness',o.Thickness,'Magnification',o.Magnification, ...
                'UsedChannels',o.UsedChannels,'StdThresh',o.StdThresh, ...
                'ConfThresh',Conf(ic),'SizeThreshUm2',SizeU(isz),'RadiusScale',Rad(ir), ...
                'OutDir',outDir,'CleanOutDir',true,'SaveOverlays',false,'SaveMasks',false, ...
                'Verbose',false); %#ok<NASGU>

            % evaluate
            R = tools.evaluate_detections(outDir, gtCsv, 0.30);
            rows(end+1,:) = { runId, Conf(ic), SizeU(isz), Rad(ir), ...
                              R.precision, R.recall, R.F1, R.tp, R.fp, R.fn, outDir }; %#ok<AGROW>
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'run','ConfThresh','SizeThreshUm2','RadiusScale','P','R','F1','TP','FP','FN','OutDir'});
T = sortrows(T, {'F1','P','R'},{'descend','descend','descend'});

% show top 10
disp(T(1:min(10,height(T)),:));

% write sweep summary
writetable(T, 'sweep_results.csv');
end
