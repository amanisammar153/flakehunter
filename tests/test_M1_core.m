classdef test_M1_core < matlab.unittest.TestCase
    methods(Test)
        function mahConfMonotonic(t)
            % Synthetic 2D RGB patch
            I = uint8( repmat(reshape(uint8(linspace(0,255,64)),1,64,1), 64,1,3) );
            p = detector.preprocess(I, struct('UsedChannels','RG','PerChannelPercentiles',[0 100]));
            % GMM: single comp centered mid-gray in used channels
            C = numel(p.usedIdx);
            mu = 0.5*ones(C,1,'single');
            Sigma = repmat(eye(C,'single')*0.05^2, 1,1,1);
            params.mu = mu; params.Sigma = Sigma; params.radius = inf(C,1,'single');
            out = detector.inferGMM(p, params, struct('Chi2Conf',true));
            % along gradient, confidence should be highest near center
            c = out.confMap(32,:); % 1D slice
            [~,imax] = max(c);
            t.verifyTrue(imax > 20 && imax < 44);
        end

        function channelSubsetParity(t)
            % Build an image where B equals mean so dropping B won't change result
            H=64; W=64;
            R = uint8(repmat(uint8(linspace(0,255,W)),H,1));
            G = R';
            B = uint8(127*ones(H,W));
            I = cat(3,R,G,B);

            % full RGB
            p1 = detector.preprocess(I, struct('UsedChannels','RGB','PerChannelPercentiles',[0 100]));
            mu = single([0.5 0.5 127/255])';
            S  = repmat(eye(3,'single')*0.1^2,1,1,1);
            params.mu = mu; params.Sigma=S; params.radius=inf(3,1,'single');
            o1 = detector.inferGMM(p1, params);

            % drop B channel; keep same mu/Sigma for RG only
            p2 = detector.preprocess(I, struct('UsedChannels','RG','PerChannelPercentiles',[0 100]));
            params2.mu = single([0.5 0.5]'); params2.Sigma = repmat(eye(2,'single')*0.1^2,1,1,1);
            params2.radius = inf(2,1,'single');
            o2 = detector.inferGMM(p2, params2);

            % Near-equality of confidences
            d = abs(single(o1.confMap) - single(o2.confMap));
            t.verifyLessThan(median(d(:)), 0.02);
        end

        function postFiltersBehave(t)
            I = uint8(zeros(128,128,3));
            p = detector.preprocess(I, struct('UsedChannels','R'));
            mu = 0.5; S = 0.05^2; params.mu = single(mu); params.Sigma = single(reshape(S,1,1,1));
            params.radius = inf(1,1,'single');
            o = detector.inferGMM(p, params);
            % fabricate a blob in conf
            cm = o.confMap; cm(40:60,40:60) = 0.9; cm(10:12,10:12) = 0.6;
            labels = uint16(ones(size(cm)));
            pp = detector.postprocess(cm, labels, struct('ConfThresh',0.7,'MinArea_um2',0,'umPerPixel',1, ...
                'MinSolidity',0.0,'KeepTopK',1,'OpenRadiusPx',1,'CloseRadiusPx',1));
            t.verifyEqual(sum(pp.mask(:)), numel(40:60)^2, 'absTol', 50); % roughly keeps the big blob
            t.verifyEqual(height(pp.stats), 1);
        end
    end
end
