function recover_ablation_tables(nSeeds)
% 作者：江苏科技大学王兴宇（Xingyu Wang）
% 功能：从已保存的原始轨迹重新计算统计量，生成论文表格与逐工况结果。
% 默认读取主程序保存的种子数量，也可显式指定需汇总的配对种子数量。
% 使用数值 CSV 写入器，避免字符串表格导出组件的兼容性问题。
rootDir=fileparts(mfilename('fullpath'));
dataDir=fullfile(rootDir,'data');
if nargin<1
    design=load(fullfile(dataDir,'ablation_design.mat'));
    if isfield(design,'nSeeds'),nSeeds=design.nSeeds;else,nSeeds=20;end
end
validateattributes(nSeeds,{'numeric'},{'scalar','real','finite','integer','positive'},mfilename,'nSeeds');
assert(mod(nSeeds,5)==0,'配对种子数量必须为正的 5 的倍数。');
names={'Coverage90_pct','Edge90_pct','TrackingRMSE90_mm','Time95_s', ...
    'Saturation90_pct','PeakDemandRatio90','MinClearance90_m','Coverage180_pct', ...
    'Seed','Q','D','F','CaseIndex'};
values=zeros(8*nSeeds,numel(names)); counter=0; zeroCount=0;
for seed=1:nSeeds
    saved=load(fullfile(dataDir,'trajectories',sprintf('seed_%02d.mat',seed)));
    p=saved.pRun;
    for variant=1:8
        tr=saved.seedTraces{variant}; win=tr.time<=p.window;
        cov=tr.coverage(find(win,1,'last'));
        edge=tr.edgeCoverage(find(win,1,'last'));
        err=1000*sqrt(mean(sum((tr.position(win,:)-tr.reference(win,:)).^2,2)));
        hit=find(tr.coverage>=95,1); elapsed=NaN;
        if ~isempty(hit),elapsed=tr.time(hit);end
        sat=100*mean(any(abs(tr.requested(win,:))>p.accelerationLimit+1e-10,2));
        ratio=max(abs(tr.requested(win,:)),[],'all')/p.accelerationLimit;
        gap=min(tr.clearance(win)); last=tr.coverage(end);
        flags=dec2bin(variant-1,3)-'0';counter=counter+1;
        values(counter,:)=[cov,edge,err,elapsed,sat,ratio,gap,last,seed,flags,saved.caseIndex];
        assert(all(isfinite(tr.position),'all'),'位置轨迹中出现非有限数值。');
        assert(all(diff(tr.coverage)>=-1e-10),'累计合格覆盖率不是单调非减序列。');
        assert(isequal(tr.position(1,:),saved.seedTraces{1}.position(1,:)),'同一种子下各组合的初始位置不一致。');
        assert(max(abs(tr.applied),[],'all')<=p.accelerationLimit+1e-12,'施加输入超出限幅范围。');
        assert(max(abs(tr.applied-max(min(tr.requested,p.accelerationLimit),-p.accelerationLimit)),[],'all')<1e-12,'施加输入与限幅后的请求输入不一致。');
        zeroCount=zeroCount+(sat==0);
    end
    fprintf('已复算并核验配对种子 %02d。\n',seed);
end
writeNumeric(fullfile(dataDir,'ablation_runs.csv'),names,values);
metrics=names(1:8);
summaryNames={'Configuration','Q','D','F','Runs','Reached95'};
for k=1:numel(metrics)
    summaryNames=[summaryNames,{[metrics{k} '_mean'],[metrics{k} '_sd']}]; %#ok<AGROW>
end
summary=zeros(8,numel(summaryNames));
byCase=zeros(40,numel(summaryNames)+1); row=0;
for variant=1:8
    flags=dec2bin(variant-1,3)-'0';
    mask=all(values(:,10:12)==flags,2);
    summary(variant,:)=summarizeGroup(values(mask,:),flags);
    for caseIndex=1:5
        row=row+1;
        byCase(row,:)=[summarizeGroup(values(mask & values(:,13)==caseIndex,:),flags),caseIndex];
    end
end
writeNumeric(fullfile(dataDir,'ablation_summary.csv'),summaryNames,summary,true);
writeNumeric(fullfile(dataDir,'ablation_by_case.csv'),[summaryNames,{'CaseIndex'}],byCase,true);
fid=fopen(fullfile(dataDir,'paired_component_contrasts.csv'),'w','n','UTF-8');
fprintf(fid,'Metric,Contrast,Mean,SD,SE\n');
metricCols=[1,2,3,5]; labels={'Q_given_DF','D_given_QF','F_given_QD','DF_interaction_given_Q'};
for col=metricCols
    matrix=reshape(values(:,col),8,nSeeds)';
    contrasts=[matrix(:,8)-matrix(:,4),matrix(:,8)-matrix(:,6),matrix(:,8)-matrix(:,7), ...
        matrix(:,8)-matrix(:,7)-matrix(:,6)+matrix(:,5)];
    for j=1:4
        fprintf(fid,'%s,%s,%.15g,%.15g,%.15g\n',names{col},labels{j},mean(contrasts(:,j)),std(contrasts(:,j)),std(contrasts(:,j))/sqrt(nSeeds));
    end
end
fclose(fid);
fid=fopen(fullfile(dataDir,'verification.txt'),'w','n','UTF-8');
fprintf(fid,'作者：江苏科技大学王兴宇（Xingyu Wang）\n已从五类工况的原始轨迹复算 %d/%d 次运行。\n限幅一致性、覆盖率单调性及配对初始条件检查通过。\n零饱和运行次数：%d。\n最小实际间距：%.9f m。\n',8*nSeeds,8*nSeeds,zeroCount,min(values(:,7)));
fclose(fid);
fid=fopen(fullfile(dataDir,'runtime.txt'),'w','n','UTF-8');fprintf(fid,'作者：江苏科技大学王兴宇（Xingyu Wang）\nMATLAB 版本：%s\n工况数量：5；配对种子数量：%d；轨迹总数：%d。\n',version,nSeeds,8*nSeeds);fclose(fid);
save(fullfile(dataDir,'ablation_verified_metrics.mat'),'values','names','summary','summaryNames','byCase','-v7');
disp(array2table(summary,'VariableNames',summaryNames));
fprintf('核验通过：%d 条轨迹；%d 次零饱和运行；最小间距 %.9f m。\n',8*nSeeds,zeroCount,min(values(:,7)));
end

function row=summarizeGroup(data,flags)
% 汇总某一组件组合的均值与样本标准差；未达到目标的时间保留 NaN。
row=[str2double(sprintf('%d%d%d',flags)),flags,size(data,1),sum(isfinite(data(:,4)))];
for col=1:8
    row=[row,mean(data(:,col),'omitnan'),std(data(:,col),'omitnan')]; %#ok<AGROW>
end
end

function writeNumeric(path,names,values,configuration)
% 第一列的组件编码保留三位格式，例如 001，不把它显示成数字 1。
if nargin<4,configuration=false;end
fid=fopen(path,'w','n','UTF-8');assert(fid>0,'无法创建 CSV 输出文件。');cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',strjoin(names,','));
for row=1:size(values,1)
    if configuration
        fprintf(fid,'%03d',values(row,1));
    else
        fprintf(fid,'%.15g',values(row,1));
    end
    for col=2:size(values,2),fprintf(fid,',%.15g',values(row,col));end
    fprintf(fid,'\n');
end
end
