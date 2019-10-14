% Problem 3.3
clc
clf
close all
clearvars

N = 255;
nfft = 1024;
fp_normx = 0.9;
fs_normx = 0.95;
f_normx = [0, fp_normx, fs_normx, 1];
ampx = [1, 1, 0, 0];
w = [1,1];
x = fir2(N,f_normx,ampx);

%% Part a
Xdft = fft(x, nfft);
f_ax = linspace(0,2,nfft);
figure()
for i = [1,2]
    ax(i) = subplot(2,1,i);
end
subplot(ax(1))
plot(f_ax,mag2db(abs(Xdft)))
xlim([0, 1])
ylim([-100, 0])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(ax(2))
plot(f_ax,(rad2deg(unwrap(angle(Xdft)))))
% plot(f_ax,(rad2deg(angle(Xdft))))
xlim([0, 1])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Phase Unwrapped (deg)')
% ylabel('Phase (deg)')
% idx_px = find(f_ax<fp_normx);
% idx_sx = intersect(find(f_ax>fs_normx),find(f_ax<1));
% dpx = max(abs(abs(Xdft(idx_px))-1));
% dsx = max(abs(abs(Xdft(idx_sx))-0));

%% Part b
% us_factor = 2;
% fp_normu = fp_normx/us_factor;
% fs_normu = fs_normx/us_factor;
% u = upsample(x,us_factor,0);
% Udft = fft(u, nfft);
% f_ax = linspace(0,2,nfft);
% figure()
% for i = [1,2]
%     ax(i) = subplot(2,1,i);
% end
% subplot(ax(1))
% plot(f_ax,mag2db(abs(Udft)))
% xlim([0, 1])
% ylim([-100, 0])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Magnitude (dB)')
% subplot(ax(2))
% plot(f_ax,(rad2deg(unwrap(angle(Udft)))))
% % plot(f_ax,(rad2deg(angle(Udft))))
% xlim([0, 1])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Phase Unwrapped (deg)')
% % ylabel('Phase (deg)')
% idx_pu = find(f_ax<fp_normu);
% idx_su = intersect(find(f_ax>fs_normu),find(f_ax<1-fs_normu));
% dpu = max(abs(abs(Udft(idx_pu))-1));
% dsu = max(abs(abs(Udft(idx_su))-0));

%% Part c
% Nf = 40;
% fp_normf = fp_normu;
% fs_normf = 2/us_factor-fp_normf;
% f_normf = [0, fp_normf, fs_normf, 1];
% ampf = [1, 1, 0, 0];
% f = firpm(Nf,f_normf,ampf);
% f = f/f(21);
% % figure()
% % stem([0:Nf],f)
% % fvtool(f)
% y = conv(u,f);
% figure()
% stem(u(256-10:257+10))
% hold on
% stem(y(276-10:277+10))
% hold off
% % fvtool(u,1,y,1)

%% Part e
% us_factor = 3;
% fp_normx2 = fp_normx/us_factor;
% fs_normx2 = fs_normx/us_factor;
% x2 = upsample(x,us_factor,0);
% 
% Nf = 80;
% fp_normf = fp_normx2;
% fs_normf = 2/us_factor-fp_normf;
% f_normf = [0, fp_normf, fs_normf, 1];
% ampf = [1, 1, 0, 0];
% f = firpm(Nf,f_normf,ampf);
% f = f/f(41);
% % figure()
% % stem([0:Nf],f)
% % fvtool(f)
% 
% z = conv(x2,f);
% Nx2 = length(x2);
% Nz = length(z);
% figure()
% stem(x2(Nx2/2-10:Nx2/2+11))
% hold on
% stem(z(Nz/2-10:Nz/2+11))
% hold off