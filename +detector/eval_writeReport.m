function REP = eval_writeReport(PR, outDir, varargin)
% DETECTOR.EVAL_WRITEREPORT  Save PR CSV + PR PNG, with an optional pick marker.
%
% REP = detector.eval_writeReport(PR, outDir, 'Prefix','report', 'PickIdx',[])
ip = inputParser; ip.CaseSensitive=false;
ip.addParameter('Prefix','report',@(s)ischar(s)||isstring(s));
ip.addParameter('PickIdx',[],@(x)isempty(x)||(isscalar(x)&&x>=1&&x<=numel(PR.thr)));
ip.parse(varargin{:});
opt = ip.Results;

if ~isfolder(outDir), mkdir(outDir); end
pref = char(opt.Prefix);

csvPath = fullfile(outDir, sprintf('%s_pr.csv', pref));
writetable(table(PR.thr(:),PR.prec(:),PR.rec(:),PR.f1(:), ...
    'VariableNames',{'threshold','precision','recall','f1'}), csvPath);

pngPath = fullfile(outDir, sprintf('%s_pr.png', pref));
figure('Visible','off'); plot(PR.rec, PR.prec, '-', 'LineWidth',1.5); grid on;
xlabel('Recall'); ylabel('Precision'); title(sprintf('PR (AP=%.3f)', PR.ap));
if ~isempty(opt.PickIdx)
    hold on; plot(PR.rec(opt.PickIdx), PR.prec(opt.PickIdx), 'o', 'MarkerSize',6, 'LineWidth',1.5);
end
saveas(gcf, pngPath); close(gcf);

REP = struct('csv',csvPath,'png',pngPath);
end
