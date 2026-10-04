function CQC_FPDFTC_ablation(nSeeds)
% 作者：江苏科技大学王兴宇（Xingyu Wang）
% 功能：运行任务空间组件消融试验，并保存可独立核验的完整轨迹。
% 三位开关从左至右依次为 Q、D、F；0 表示未启用，1 表示启用。
% Q：利用已观测的表面质量选择观测位置与方向。
% D：延迟激活的非线性跟踪反馈及基于误差的参考轨迹进度调节。
% F：采用归一化高斯基函数在线估计扰动。
% 各配置共用同一个简化笛卡尔伺服对象；这不是完整关节力矩模型。
% 调用示例：CQC_FPDFTC_ablation(20)，对应五类工况、160 条轨迹。
if nargin < 1, nSeeds = 20; end
validateattributes(nSeeds,{'numeric'},{'scalar','real','finite','integer','positive'},mfilename,'nSeeds');
assert(mod(nSeeds,5)==0,'配对种子数量必须为正的 5 的倍数，以保证五类工况数量相同。');
author='江苏科技大学王兴宇（Xingyu Wang）';
rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir, 'data');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
traceDir = fullfile(dataDir, 'trajectories');
if ~exist(traceDir, 'dir'), mkdir(traceDir); end
p = parameters();
flags = dec2bin(0:7, 3) - '0';
for seed = 1:nSeeds
    seedTraces = cell(1,8);
    pRun=p;
    caseIndex=mod(seed-1,5)+1;
    if ismember(caseIndex,[2,5])
        pRun.depth=1.18*p.depth; pRun.twist=1.10*p.twist;
    end
    if ismember(caseIndex,[3,5]), pRun.initialOffset=[0.025,-0.018,0.020]; end
    if ismember(caseIndex,[4,5])
        pRun.obstacles=[p.obstacles;0.12,0.42,0.30,0.17];
    end
    plate=makePlate(pRun);
    rng(7300 + seed, 'twister');
    noise = randn(round(p.duration/p.dt)+1, 3);
    disturbancePhase = 2*pi*rand(1, 3);
    payload = 1.1 + 0.15*rand;
    for k = 1:8
        [~, trace] = runCase(flags(k,:), pRun, plate, noise, disturbancePhase, payload);
        seedTraces{k} = trace;
    end
    save(fullfile(traceDir,sprintf('seed_%02d.mat',seed)), 'seedTraces', 'seed', 'payload', 'disturbancePhase', 'pRun', 'caseIndex', '-v7.3');
    fprintf('已完成配对种子 %d/%d。\n', seed, nSeeds);
end
save(fullfile(dataDir, 'ablation_design.mat'), 'p', 'flags', 'nSeeds', 'author', '-v7');
recover_ablation_tables(nSeeds);
end

function p = parameters()
% 所有八种组件组合使用同一组参数；各工况的改动统一作用于该种子内的所有组合。
p.length = 2.80; p.width = 1.02; p.depth = 0.42; p.twist = 0.055;
p.Nu = 72; p.Nv = 46; p.standoff = 0.46;
p.qualityThreshold = 0.52; p.footprint = 0.245; p.rangeSigma = 0.080;
p.maxIncidence = 52; p.safeDistance = 0.14;
p.dt = 0.02; p.frameStep = 5; p.duration = 180; p.window = 90;
p.speed = 0.36; p.accelerationLimit = 1.2; p.referenceAccel = 0.65;
p.kp = 49; p.kd = 14; p.damping = 0.35;
p.sensorNoise = 0.00012; p.bound = 0.012;
p.initialOffset=[0.006,-0.005,0.008];
p.learningRate = 3.0; p.leakage = 0.025; p.weightBound = 0.7;
p.obstacles = [-0.78,0.58,0.30,0.15;0.86,-0.56,0.36,0.16;1.36,0.20,0.42,0.13];
end

function plate = makePlate(p)
u = linspace(-1,1,p.Nu); v = linspace(-1,1,p.Nv);
[U,V] = meshgrid(u,v);
X = p.length*U/2; Y = p.width*V/2;
Z = -0.36+p.depth*U.^2+p.twist*sin(pi*V).*cos(pi*U/2);
zx = (2*p.depth*U - p.twist*pi/2*sin(pi*V).*sin(pi*U/2))*2/p.length;
zy = p.twist*pi*cos(pi*V).*cos(pi*U/2)*2/p.width;
N = [-zx(:),-zy(:),ones(numel(U),1)];
N = N ./ vecnorm(N,2,2);
plate.points = [X(:),Y(:),Z(:)]; plate.normals = N;
plate.edge = max(abs(U(:)),abs(V(:))) >= 0.74;
plate.U = U; plate.V = V;
uv = zeros(18*17,2);
for j=1:18
    vv = linspace(-0.94,0.94,17);
    if mod(j,2)==0, vv=fliplr(vv); end
    uv((j-1)*17+(1:17),:) = [repmat(-0.96+(j-1)*1.92/17,17,1),vv'];
end
plate.target = zeros(size(uv,1),3);
plate.targetNormal = zeros(size(uv,1),3);
for j=1:size(uv,1)
    a=uv(j,1); b=uv(j,2);
    point=[p.length*a/2,p.width*b/2,-0.36+p.depth*a^2+p.twist*sin(pi*b)*cos(pi*a/2)];
    normal=[-(2*p.depth*a-p.twist*pi/2*sin(pi*b)*sin(pi*a/2))*2/p.length, ...
        -p.twist*pi*cos(pi*b)*cos(pi*a/2)*2/p.width,1];
    normal=normal/norm(normal);
    sensor=point+p.standoff*normal;
    for pass=1:6
        for o=1:size(p.obstacles,1)
            dvec=sensor-p.obstacles(o,1:3); minD=p.obstacles(o,4)+p.safeDistance+0.060;
            if norm(dvec)<minD, sensor=p.obstacles(o,1:3)+minD*dvec/norm(dvec); end
        end
    end
    plate.target(j,:)=sensor; plate.targetNormal(j,:)=normal;
end
% 候选视点质量只由几何和观测条件计算，不按组件组合名称额外加分或扣分。
plate.viewScores=zeros(numel(U),size(plate.target,1));
for j=1:size(plate.target,1)
    plate.viewScores(:,j)=viewQuality(plate.target(j,:),plate.targetNormal(j,:),plate,p);
end
end

function [out,tr] = runCase(flags,p,plate,noise,phase,payload)
n=round(p.duration/p.dt)+1; time=(0:n-1)'*p.dt;
position=zeros(n,3); reference=position; requested=position; applied=position;
quality=zeros(size(plate.points,1),1); coverage=zeros(n,1); edgeCoverage=coverage;
errorHistory=coverage; clearance=coverage; selected=zeros(n,1);
x=plate.target(1,:)+p.initialOffset; velocity=[0,0,0];
xRef=plate.target(1,:); vRef=[0,0,0]; oldVRef=vRef;
axisNow=plate.targetNormal(1,:); current=1;
centers=[-1,-1;-1,0;-1,1;0,-1;0,0;0,1;1,-1;1,0;1,1];
weights=zeros(9,3); previousBasis=ones(9,1)/9;
previousVelocity=velocity; previousApplied=[0,0,0];
sampledQuality=0; sampledEdge=0;
for it=1:n
    t=time(it);
    if flags(3) && it>1
        observed=(velocity-previousVelocity)/p.dt-previousApplied+p.damping*previousVelocity;
        prediction=previousBasis'*weights;
        innovation=observed-prediction;
        weights=weights+p.dt*(p.learningRate*previousBasis*innovation/(0.05+sum(previousBasis.^2))-p.leakage*weights);
        weights=max(min(weights,p.weightBound),-p.weightBound);
    end
    if norm(xRef-plate.target(current,:))<0.025
        if flags(1)
            deficiency=max(0,p.qualityThreshold-quality);
            gain=(deficiency' * plate.viewScores)';
            travel=vecnorm(plate.target-xRef,2,2);
            gain(travel<0.08)=0;
            utility=gain./(0.35+travel);
            if max(utility)>1e-8
                [~,current]=max(utility);
            else
                current=mod(current,size(plate.target,1))+1;
            end
        else
            current=mod(current,size(plate.target,1))+1;
        end
    end
    direction=plate.target(current,:)-xRef;
    desiredV=p.speed*direction/max(norm(direction),0.02);
    e=xRef-x;
    activation=smoothStep(min(max((t-1)/2,0),1));
    if flags(2)
        governor=max(0.15,min(1,p.bound/(norm(e)+0.002)));
        desiredV=governor*desiredV;
    end
    dv=desiredV-vRef;
    vRef=vRef+dv*min(1,p.referenceAccel*p.dt/max(norm(dv),eps));
    % 两个目标位置安全并不代表它们之间的连线安全。
    % 因此所有组合共用参考轨迹可行性过滤器，先约束参考运动本身。
    for pass=1:3
        for o=1:size(p.obstacles,1)
            radial=xRef-p.obstacles(o,1:3); normal=radial/norm(radial);
            gap=norm(radial)-p.obstacles(o,4)-p.safeDistance-0.045;
            lowerRate=-3*gap;
            if dot(vRef,normal)<lowerRate
                vRef=vRef+(lowerRate-dot(vRef,normal))*normal;
            end
        end
    end
    xRef=xRef+p.dt*vRef;
    aRef=(vRef-oldVRef)/p.dt;
    measuredX=x+p.sensorNoise*noise(it,:);
    e=xRef-measuredX; ev=vRef-velocity;
    basis=exp(-sum((centers-[x(1)/1.5,x(2)/0.6]).^2,2)/(2*0.7^2));
    basis=basis/sum(basis);
    compensation=[0,0,0];
    if flags(3), compensation=basis'*weights; end
    u=aRef+p.kp*e+p.kd*ev+p.damping*velocity-compensation;
    if flags(2)
        scaled=e/p.bound;
        u=u+activation*(0.10*sign(scaled).*abs(scaled).^0.6+0.04*sign(scaled).*abs(scaled).^1.4);
    end
    % 所有组合使用相同的夹具安全约束，再统一进行实际输入限幅。
    for o=1:size(p.obstacles,1)
        radial=x-p.obstacles(o,1:3); radius=norm(radial);
        normal=radial/radius; gap=radius-p.obstacles(o,4)-p.safeDistance;
        bound=-7*dot(velocity,normal)-16*gap;
        if dot(u,normal)<bound, u=u+(bound-dot(u,normal))*normal; end
    end
    uApplied=max(min(u,p.accelerationLimit),-p.accelerationLimit);
    position(it,:)=x; reference(it,:)=xRef; requested(it,:)=u; applied(it,:)=uApplied;
    errorHistory(it)=norm(xRef-x);
    clearance(it)=min(vecnorm(p.obstacles(:,1:3)-x,2,2)-p.obstacles(:,4));
    selected(it)=current;
    if flags(1), targetAxis=plate.targetNormal(current,:); else, targetAxis=[0,0,1]; end
    axisNow=axisNow+min(1,4*p.dt)*(targetAxis-axisNow); axisNow=axisNow/norm(axisNow);
    if mod(it-1,p.frameStep)==0
        q=viewQuality(x,axisNow,plate,p);
        % 有限曝光造成的运动模糊取决于实际速度；各组合使用同一个衰减模型。
        q=q*exp(-0.18*norm(velocity)^2);
        quality=max(quality,q);
        sampledQuality=100*mean(quality>=p.qualityThreshold);
        sampledEdge=100*mean(quality(plate.edge)>=p.qualityThreshold);
    end
    coverage(it)=sampledQuality; edgeCoverage(it)=sampledEdge;
    disturbance=0.16*sin(0.11*t+phase)+[0.14,-0.10,0.12] ...
        +0.22*exp(-((t-42)/5)^2)*[1,-0.6,0.8];
    previousVelocity=velocity; previousApplied=uApplied; previousBasis=basis;
    acceleration=uApplied/payload-p.damping*velocity+disturbance;
    velocity=velocity+p.dt*acceleration;
    x=x+p.dt*velocity;
    oldVRef=vRef;
end
window=time<=p.window;
out.Coverage90_pct=coverage(find(window,1,'last'));
out.Edge90_pct=edgeCoverage(find(window,1,'last'));
out.TrackingRMSE90_mm=1000*sqrt(mean(errorHistory(window).^2));
hit=find(coverage>=95,1);
if isempty(hit), out.Time95_s=NaN; else, out.Time95_s=time(hit); end
out.Saturation90_pct=100*mean(any(abs(requested(window,:))>p.accelerationLimit+1e-10,2));
out.PeakDemandRatio90=max(abs(requested(window,:)),[],'all')/p.accelerationLimit;
out.MinClearance90_m=min(clearance(window));
out.Coverage180_pct=coverage(end);
out.Seed=0; out.Q=0; out.D=0; out.F=0;
tr=struct('time',time,'position',position,'reference',reference,'requested',requested, ...
    'applied',applied,'coverage',coverage,'edgeCoverage',edgeCoverage,'error',errorHistory, ...
    'clearance',clearance,'selectedView',selected,'qualityField',quality);
assert(all(abs(applied)<=p.accelerationLimit+1e-12,'all'),'施加输入超出声明的限幅范围。');
assert(all(diff(coverage)>=-1e-12),'累计合格覆盖率出现了非物理的下降。');
assert(all(isfinite(position),'all'),'位置轨迹中出现非有限数值。');
end

function q = viewQuality(position,axis,plate,p)
delta=plate.points-position;
distance=vecnorm(delta,2,2);
cosInc=max(0,plate.normals*axis');
axial=-delta*axis';
lateral2=max(0,sum(delta.^2,2)-axial.^2);
q=exp(-lateral2/(2*p.footprint^2)).*cosInc.^1.7.*exp(-(distance-p.standoff).^2/(2*p.rangeSigma^2));
q(cosInc<cosd(p.maxIncidence) | axial<=0)=0;
% 对观测射线段与球形排除区域进行相交检测，该检测与组件名称无关。
for j=1:size(p.obstacles,1)
    oc=p.obstacles(j,1:3)-position;
    parameter=(delta*oc')./max(sum(delta.^2,2),eps);
    closest=position+max(0,min(1,parameter)).*delta;
    hit=sum((closest-p.obstacles(j,1:3)).^2,2)<p.obstacles(j,4)^2;
    q(hit & parameter>0 & parameter<1)=0;
end
q=min(1,max(0,q));
end

function y=smoothStep(x)
y=3*x^2-2*x^3;
end
