clear,clc;
%% 参数定义
c = 3e8;
% 波形
wave.f0 = 20e9;
wave.B = 10e6;
wave.prf = 40e3;
wave.K = 32;     % 脉冲个数
wave.duty_cycle = 0.2;
wave.prt = wave.duty_cycle/wave.prf;
wave.pri = 1/wave.prf;
wave.lambda = c/wave.f0;
wave.mu = wave.B/wave.prt;   %调频斜率
% 采样率
fs = 80e6;
% 平台位置速度
plat.position = [0;0;1e3];
plat.v = [0;100;0];
% 探测距离
min_range = c/2*wave.prt;
max_range = c/2*wave.pri;
t = (0:wave.pri*fs-1)/fs;

%% 发射阵列
Trans.M = 16;
Trans.N = 16;
Trans.Pt = 5e3;
Trans.Gt_db = 0;
Trans.d = wave.lambda/2;
Trans.Gt = 10^(Trans.Gt_db/10);
Trans.angle = [[49.9;-26],[29;-25],[10;-24],[39;-23]];
% Trans.angle = [60;-50];
Trans.w = zeros(Trans.M*Trans.N,1);

for i = 1:size(Trans.angle,2)
    Trans.w = Trans.w + array_space(Trans.M,Trans.N,Trans.d,wave.lambda,Trans.angle(1,i),Trans.angle(2,i));
end
Trans.w = conj(Trans.w)/size(Trans.angle,2);

% ant_grid_n = 181;
% angel_grid = linspace(-90,90,ant_grid_n);
% ant_gra = zeros(ant_grid_n,ant_grid_n);
% for i = 1:ant_grid_n
%     for j = 1:ant_grid_n
%         ant_gra(i,j) = abs(Trans.w.'*array_space(Trans.M,Trans.N,Trans.d,wave.lambda,angel_grid(i),angel_grid(j)));
%     end
% end
% mesh(angel_grid,angel_grid,ant_gra);

%% 接收阵列
Rece.M = Trans.M;
Rece.N = Trans.N;
Rece.d = Trans.d;
Rece.Gr = Trans.Gt;

%% 噪声
noise.F = 1;
noise.NF = 10*log10(noise.F);
noise.dbm = -174 + 10*log10(wave.B) + noise.NF;
noise.w = 10^((noise.dbm-30)/10);

%% 目标回波 列表示第几个目标
tgt.pos = [[1262.9,1505,10].',[1797.9,1038,10].',[2153.1,379.6,10].',[1758.6,1475.6,10].'];
tgt.v = [[0,-20,0].',[0,-20,0].',[0,-20,0].',[0,-20,0].'];
tgt.rcs = [0.01,0.01,0.01,0.01];
tgt.num = length(tgt.pos);

for i = 1:tgt.num
    tgt.dis(i) = norm(plat.position-tgt.pos(:,i));
end

for i = 1:tgt.num
    tgt.angle(:,i) = compute_angles(plat.position, tgt.pos(:,i));
end

for i = 1:tgt.num
    tgt.amp(i) = sqrt(Trans.Pt*Trans.Gt*Rece.Gr*wave.lambda^2*tgt.rcs(i)/((4*pi)^3*tgt.dis(i)^4))*...
        Trans.w.'*array_space(Trans.M,Trans.N,Trans.d,wave.lambda,tgt.angle(1,i),tgt.angle(2,i));
end

for i = 1:tgt.num
    tgt.tao0(i) = 2*tgt.dis(i)/c;
end

for i = 1:tgt.num
    tgt.vd(i) = (tgt.v(:,i) - plat.v).'*(tgt.pos(:,i)-plat.position)/norm(tgt.pos(:,i)-plat.position);
    tgt.fd(i) = -2/wave.lambda*tgt.vd(i);
end

for i = 1:tgt.num
    tgt.fu_z(:,i) = Rece.d/wave.lambda*sind(tgt.angle(2,i));
    tgt.fu_y(:,i) = Rece.d/wave.lambda*cosd(tgt.angle(2,i))*sind(tgt.angle(1,i));
end

% % 生成目标回波信号 
tgt.data = zeros(Rece.N*Rece.M*wave.K,wave.pri*fs);
for i = 1:tgt.num
        tgt.a_u = kron(exp(-1j*2*pi*(0:Rece.N-1)*tgt.fu_y(i)).',exp(1j*2*pi*(0:Rece.M-1)*tgt.fu_z(i)).');
        tgt.a_d = exp(-1j*2*pi*(0:wave.K-1)*tgt.fd(i)*wave.pri).';
        tgt.a_st = kron(tgt.a_u,tgt.a_d);
        tgt.data = tgt.data + tgt.a_st*tgt.amp(i)*rectpuls((t-tgt.tao0(i))/wave.prt,1).*...
            exp(1j*pi*wave.mu*(t-tgt.tao0(i)).^2).*exp(1j*2*pi*tgt.fd(i)*t);
end

%% 回波合成
rng(150);
echo.data = sqrt(noise.w)*(randn(size(tgt.data)) + 1j*randn(size(tgt.data))); % 噪声
echo.data = echo.data + tgt.data;
echo.dis_gate = (0:size(echo.data,2)-1)*c/2/fs;
echo.v_gate = (-wave.K/2:wave.K/2-1)*wave.lambda/2/wave.K*wave.prf;
% echo.v_gate = linspace(-wave.prf/2,wave.prf/2,wave.K);
% plot(real(echo.data(1,:)));

%% 脉冲压缩
pc.pri_num = wave.pri*fs;
pc.prt_num = wave.prt*fs;
pc.nfft = 2^nextpow2(pc.pri_num + pc.prt_num - 1);
pc.delay = floor(pc.prt_num/2);
pc.t = (-wave.prt*fs/2:wave.prt*fs/2-1)/fs;
pc.flter_index = rectpuls(-pc.t/wave.prt,1).*exp(-1j*pi*wave.mu*pc.t.^2);
pc.flter_w = fft(pc.flter_index,pc.nfft);

pc.out = fft(echo.data,pc.nfft,2);
for i = 1:size(pc.out,1)
    pc.out(i,:) = pc.out(i,:).*pc.flter_w;
end
pc.out = ifft(pc.out,pc.nfft,2);
pc.out = pc.out(:,pc.delay+1:pc.delay+pc.pri_num);
% plot(echo.dis_gate,abs(pc.out(1,:)));

%% 波束形成
beamform.angle(1) = Trans.angle(1,1);
beamform.angle(2) = Trans.angle(2,1);
beamform.w = conj(array_space(Trans.M,Trans.N,Trans.d,wave.lambda,beamform.angle(1),beamform.angle(2)));
for i = 1:wave.K
    bufer_one = pc.out(i:wave.K:end,:);
    beamform.out(i,:) = beamform.w.'*bufer_one;
end
% plot(real(beamform.out(1,:)));

%% MTD
mtd.out = ifftshift(fft(bsxfun(@times, beamform.out, hanning(wave.K))),1);
% mesh(echo.dis_gate,echo.v_gate,10*log10((abs(mtd.out).^2)));
% histogram(abs(mtd.out(1,:)).^2);

%% CFAR
[cfar_out,cfar_ref] = cfar_ca(mtd.out,1e-6,20,3);
% mesh(echo.dis_gate,echo.v_gate,cfar_out);
% plot(abs(mtd.out(8,:).^2));
% hold on;
% plot(cfar_ref(8,:));


%% 聚类
cluster_tgt.index = find(any(cfar_out));
cluster_tgt.num = length(cluster_tgt.index);
cluster_tgt.num_line = ones(length(cluster_tgt.index),1); % 每个距离单元目标数量

%% doa
doa.angle = zeros(2,cluster_tgt.num);

for i = 1:cluster_tgt.num
    doa.input = reshape(pc.out(:,cluster_tgt.index(i)),wave.K,Trans.M*Trans.N).';
    doa.rx = 1/wave.K*(doa.input*doa.input');
    [doa.eig_vec,doa.eig_val_d] = eig(doa.rx);
    doa.eig_val = diag(doa.eig_val_d);
    [doa.eig_val,doa.sort_index] = sort(doa.eig_val); %特征值从小到大排序
    doa.eig_vec = doa.eig_vec(:,doa.sort_index);
    doa.ns = doa.eig_vec(:,1:end-cluster_tgt.num_line(i)); %生成噪声空间
    % 谱搜索
    doa.grid_n = 181;
    doa.grid = linspace(-90,90,doa.grid_n); %搜索间隔
    doa.music = zeros(doa.grid_n,doa.grid_n);
    for j = 1:doa.grid_n
        for k = 1:doa.grid_n
            doa.a_s = array_space(Trans.M,Trans.N,Trans.d,wave.lambda,doa.grid(j),doa.grid(k));
            doa.music(j,k) = 1/((doa.ns'*doa.a_s)'*(doa.ns'*doa.a_s));
        end
    end
    % mesh(doa.grid,doa.grid,doa.music);
    [doa.max_row,doa.max_col] = find(doa.music == max(max(doa.music)));
    doa.angle(1,i) = doa.grid(doa.max_row);
    doa.angle(2,i) = doa.grid(doa.max_col);
end

%% 关联数据
dect_tgt.dis = echo.dis_gate(cluster_tgt.index)
dect_tgt.angel = doa.angle