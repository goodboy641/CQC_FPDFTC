function verify_ablation()
% 作者：江苏科技大学王兴宇（Xingyu Wang）
% 功能：独立复算原始轨迹中的关键指标，并与未取整的 CSV 数值逐项比较。
% 同时核验组件顺序、工况编号、限幅一致性及配对初始条件。
rootDir=fileparts(mfilename('fullpath'));
dataDir=fullfile(rootDir,'data');
raw=readtable(fullfile(dataDir,'ablation_runs.csv'));
design=load(fullfile(dataDir,'ablation_design.mat')); p=design.p;
if isfield(design,'nSeeds'),nSeeds=design.nSeeds;else,nSeeds=20;end
assert(height(raw)==8*nSeeds,'CSV 行数与设计的运行次数不一致。');
assert(isequal(unique(raw.Seed),(1:nSeeds)'),'CSV 中的配对种子编号不完整。');
zeroRuns=0; lowestClearance=inf;
for seed=1:nSeeds
    saved=load(fullfile(dataDir,'trajectories',sprintf('seed_%02d.mat',seed)));
    group=raw(raw.Seed==seed,:);
    assert(height(group)==8,'每个配对种子必须包含八种组件组合。');
    assert(isequal([group.Q,group.D,group.F],dec2bin(0:7,3)-'0'),'CSV 的组件顺序必须为 Q、D、F，组合从 000 到 111。');
    assert(all(group.CaseIndex==saved.caseIndex),'CSV 工况编号与轨迹文件不一致。');
    for k=1:8
        tr=saved.seedTraces{k}; win=tr.time<=p.window;
        sat=100*mean(any(abs(tr.requested(win,:))>p.accelerationLimit+1e-10,2));
        err=1000*sqrt(mean(sum((tr.reference(win,:)-tr.position(win,:)).^2,2)));
        assert(abs(sat-group.Saturation90_pct(k))<1e-10,'饱和率复算结果与 CSV 不一致。');
        assert(abs(err-group.TrackingRMSE90_mm(k))<1e-8,'跟踪误差复算结果与 CSV 不一致。');
        assert(max(abs(tr.applied),[],'all')<=p.accelerationLimit+1e-12,'施加输入超出限幅范围。');
        assert(max(abs(tr.applied-max(min(tr.requested,p.accelerationLimit),-p.accelerationLimit)),[],'all')<1e-12,'施加输入与请求输入的限幅结果不一致。');
        assert(abs(tr.coverage(find(win,1,'last'))-group.Coverage90_pct(k))<1e-10,'90 秒覆盖率复算结果不一致。');
        assert(isequal(tr.position(1,:),saved.seedTraces{1}.position(1,:)),'配对组合的初始位置不一致。');
        hit=find(tr.coverage>=95,1);
        if isempty(hit)
            assert(isnan(group.Time95_s(k)),'未达到 95% 覆盖率的运行必须保留 NaN。');
        else
            assert(abs(tr.time(hit)-group.Time95_s(k))<1e-10,'首次达到 95% 覆盖率的时间不一致。');
        end
        zeroRuns=zeroRuns+(sat==0);
        lowestClearance=min(lowestClearance,min(tr.clearance(win)));
    end
end
fprintf('核验通过：%d 条轨迹；限幅、跟踪误差、覆盖率、未达标记录及配对初始条件一致。\n',8*nSeeds);
fprintf('零饱和运行次数：%d / %d；最小观测间距 %.9f m。\n',zeroRuns,8*nSeeds,lowestClearance);
fid=fopen(fullfile(dataDir,'verification.txt'),'w','n','UTF-8');
fprintf(fid,'作者：江苏科技大学王兴宇（Xingyu Wang）\n%d/%d 条轨迹的指标复算全部通过。\nQ、D、F 顺序及工况编号检查通过。\n零饱和运行次数：%d。\n最小观测间距：%.9f m。\n',8*nSeeds,8*nSeeds,zeroRuns,lowestClearance);
fclose(fid);
end
