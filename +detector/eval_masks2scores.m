classdef test_M3_evalPR < matlab.unittest.TestCase
    methods(Test)
        function pr_behaves_monotonic(t)
            H=80; W=120;
            [X,Y] = meshgrid(linspace(-1,1,W),linspace(-1,1,H));
            r = sqrt(X.^2+Y.^2);
            conf = single(1 - min(1,max(0,r)));   % high in center, low at edges
            gt   = r <= 0.4;                      % disk ground truth

            [sp,sn] = detector.eval_masks2scores(conf, gt);
            PR = detector.eval_computePR(sp, sn, 'Thresh', linspace(0,1,51));

            % sanity: recall drops (or stays) as threshold increases
            t.verifyTrue(all(diff(PR.rec) <= 1e-9));
            % precision usually increases as threshold increases (not strictly guaranteed)
            t.verifyGreaterThanOrEqual(PR.prec(end), PR.prec(1));
            % AP should be reasonable (>0.5 for this clean case)
            t.verifyGreaterThan(PR.ap, 0.5);
        end
    end
end
