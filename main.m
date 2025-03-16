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
Trans.N = 1;
Trans.M = 1;
Trans.Pt = 50e3;
Trans.Gt_db = 20;
Trans.Gt = 10^(Trans.Gt_db/10);

%% 接收阵列
Rece.N = 1;
Rece.M = 1;
Rece.d = wave.lambda/2;
Rece.Gr = Trans.Gt;

%% 噪声
noise.F = 1;
noise.NF = 10*log10(noise.F);
noise.dbm = -174 + 10*log10(wave.B) + noise.NF;
noise.w = 10^((noise.dbm-30)/10);

%% 目标回波 列表示第几个目标
tgt.pos = [[1262.9,1505,10].',[1797.9,1038,10].',[2153.1,379.6,10].',[1758.6,1475.6,10].'];
tgt.v = [[0,-20,0].',[0,-20,0].',[0,-20,0].',[0,-20,0].'];
tgt.rcs = [1,1,1,1];
tgt.num = length(tgt.pos);

for i = 1:tgt.num
    tgt.dis(i) = norm(plat.position-tgt.pos(:,i));
end

for i = 1:tgt.num
    tgt.amp(i) = sqrt(Trans.Pt*Trans.Gt*Rece.Gr*wave.lambda^2*tgt.rcs(i)/((4*pi)^3*tgt.dis(i)^4));
end

for i = 1:tgt.num
    tgt.tao0(i) = 2*tgt.dis(i)/c;
end

for i = 1:tgt.num
    tgt.angle(:,i) = compute_angles(plat.position, tgt.pos(:,i));
end

for i = 1:tgt.num
    tgt.fu(:,i) = Rece.d/wave.lambda*cosd(tgt.angle(1,i))*cosd(tgt.angle(2,i));
end

for i = 1:tgt.num
    tgt.fd(i) = 2/wave.lambda*(plat.v - tgt.v(:,i)).'*(tgt.pos(:,i)-plat.position)/norm(tgt.pos(:,i)-plat.position);
end

% 生成目标回波信号 
tgt.data = zeros(Rece.N*Rece.M*wave.K,wave.pri*fs);
for i = 1:tgt.num
        tgt.a_u = exp(-1j*2*pi*(0:Rece.N-1)*tgt.fu(i)).';
        tgt.a_d = exp(-1j*2*pi*(0:wave.K-1)*tgt.fd(i)*wave.pri).';
        tgt.a_st = kron(tgt.a_u,tgt.a_d);
        tgt.data = tgt.data + tgt.a_st*tgt.amp(i)*rectpuls((t-tgt.tao0(i))/wave.prt,1).*...
            exp(1j*pi*wave.mu*(t-tgt.tao0(i)).^2).*exp(1j*2*pi*tgt.fd(i)*t);
end

%% 杂波
% range_resolution = c/2/wave.B;
% degree_resolution = 3;
% clutter.start = 0;
% clutter.end = sqrt(max_range^2-plat.position(3)^2);
% clutter.sigma0_dbm = -20;
% clutter.sigma0 = 10^(clutter.sigma0_dbm/10);
% clutter.dis_num = floor((clutter.end - clutter.start)/range_resolution);
% clutter.deg_num = 180/degree_resolution;
% clutter.center_pos = zeros(clutter.deg_num,clutter.dis_num,3);
% clutter.center_dis = zeros(clutter.deg_num,clutter.dis_num);
% clutter.center_theta = zeros(clutter.deg_num,clutter.dis_num);
% clutter.center_phi = zeros(clutter.deg_num,clutter.dis_num);
% clutter.center_size = zeros(clutter.deg_num,clutter.dis_num);
% clutter.rcs = zeros(clutter.deg_num,clutter.dis_num);
% clutter.Kr = zeros(clutter.deg_num,clutter.dis_num);
% clutter.fu = zeros(clutter.deg_num,clutter.dis_num);
% clutter.fd = zeros(clutter.deg_num,clutter.dis_num);
% clutter.tao0 = zeros(clutter.deg_num,clutter.dis_num);
% clutter.data = zeros(Rece.N*Rece.M*wave.K,wave.pri*fs);
% 
% for i = 1:clutter.deg_num
%     for j = 1:clutter.dis_num
%         clutter.center_pos(i,j,1) = range_resolution*(j-0.5)*sind(-90+(i-0.5)*degree_resolution); 
%         clutter.center_pos(i,j,2) = range_resolution*(j-0.5)*cosd(-90+(i-0.5)*degree_resolution); 
%     end
% end
% 
% for i = 1:clutter.deg_num
%     for j = 1:clutter.dis_num
%         clutter.center_dis(i,j) = norm(plat.position-squeeze(clutter.center_pos(i,j,:)));
%     end
% end
% clutter.center_theta = repmat((((1:clutter.deg_num)-0.5)*degree_resolution - 90).',1,clutter.dis_num);
% 
% for j = 1:clutter.dis_num
%     clutter.center_phi(1,j) = atan2d(norm(plat.position(3)-squeeze(clutter.center_pos(1,j,3)))...
%         ,norm(plat.position(1:2)-squeeze(clutter.center_pos(1,j,1:2))));
% end
% clutter.center_phi = repmat(clutter.center_phi(1,:),clutter.deg_num,1);
% 
% for i = 1:clutter.dis_num
%     clutter.center_size(1,i) = degree_resolution/360*(pi*(i*range_resolution)^2-pi*((i-1)*range_resolution)^2);
% end
% clutter.center_size = repmat(clutter.center_size(1,:),clutter.deg_num,1);
% 
% clutter.rcs = clutter.sigma0*clutter.center_size;
% 
% for i = 1:clutter.dis_num
%     clutter.Kr(1,i) = sqrt(Trans.Pt*Trans.Gt*Rece.Gr*wave.lambda^2*clutter.rcs(1,i)...
%         /((4*pi)^3*clutter.center_dis(1,i)^4));
% end
% clutter.Kr = repmat(clutter.Kr(1,:),clutter.deg_num,1);
% 
% clutter.fu = Rece.d/wave.lambda*cosd(clutter.center_theta).*cosd(clutter.center_phi);
% 
% clutter.fd = 2*plat.v(2)/wave.lambda*cosd(clutter.center_theta).*cosd(clutter.center_phi);
% 
% for i = 1:clutter.dis_num
%     clutter.tao0(1,i) = 2*clutter.center_dis(1,i)/c;
% end
% clutter.tao0 = repmat(clutter.tao0(1,:),clutter.deg_num,1);
% 
% % 生成杂波
% for i = 1:clutter.deg_num
%     for j = 1:clutter.dis_num
%         a_u = exp(-1j*2*pi*(0:Rece.N-1)*clutter.fu(i,j)).';
%         a_d = exp(-1j*2*pi*(0:wave.K-1)*clutter.fd(i,j)*wave.pri).';
%         a_st = kron(a_u,a_d);
%         clutter.data = clutter.data + a_st*clutter.Kr(i,j)*rectpuls((t-clutter.tao0(i,j))/wave.prt,1).*...
%             exp(1j*pi*wave.mu*(t-clutter.tao0(i,j)).^2).*exp(1j*2*pi*clutter.fd(i,j)*t);
%     end
% end

%% 回波合成
echo.data = sqrt(noise.w)*(randn(size(tgt.data)) + 1j*randn(size(tgt.data))); % 噪声
% echo.real_noise_power = 1/size(echo.data,2)/2*sum(diag(echo.data*echo.data'))/size(echo.data,1);
echo.data = echo.data + tgt.data;
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
% plot((0:length(pc.out)-1)*c/2/fs,abs(pc.out(1,:)));

%% MTD
% mesh(10*log10(abs(ifftshift(fft(bsxfun(@times, pc.out, hanning(32))),1))));
mtd.out = ifftshift(fft(bsxfun(@times, pc.out, hanning(wave.K))),1);
% mesh(10*log10((abs(mtd.out).^2)));
% histogram(abs(mtd.out(1,:)).^2);

%% CFAR
[cfar_out,cfar_ref] = cfar_ca(mtd.out,1e-6,20,3);
mesh(cfar_out);
% plot(abs(mtd.out(8,:).^2));
% hold on;
% plot(cfar_ref(8,:));