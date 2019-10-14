% Problem 3.4
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
x1 = downsample(x,2,0);
% X1dft = fft(x1, nfft);
% f_ax = linspace(0,2,nfft);
% figure()
% for i = [1,2]
%     ax(i) = subplot(2,1,i);
% end
% subplot(ax(1))
% plot(f_ax,mag2db(abs(X1dft)))
% xlim([0, 1])
% ylim([-100, 0])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Magnitude (dB)')
% subplot(ax(2))
% plot(f_ax,(rad2deg(unwrap(angle(X1dft)))))
% % plot(f_ax,(rad2deg(angle(X1dft))))
% xlim([0, 1])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Phase Unwrapped (deg)')
% ylabel('Phase (deg)')
% idx_px = find(f_ax<fp_normx);
% idx_sx = intersect(find(f_ax>fs_normx),find(f_ax<1));
% dpx = max(abs(abs(X1dft(idx_px))-1));
% dsx = max(abs(abs(X1dft(idx_sx))-0));

%% Part b
Nh1 = 80;
ds_factor = 2;
fc_normh1 = 1/ds_factor;
h1 = fir1(Nh1,fc_normh1);
h1 = h1/h1(Nh1/2+1);
% figure()
% stem(h1)
% fvtool(h1)
y1 = conv(x,h1);
% fvtool(y1)

%% Part c
v1 = downsample(y1,2,0);
V1dft = fft(v1,nfft);
f_ax = linspace(0,2,nfft);
figure()
for i = [1,2]
    ax(i) = subplot(2,1,i);
end
subplot(ax(1))
plot(f_ax,mag2db(abs(V1dft)))
xlim([0, 1])
ylim([-100, 0])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(ax(2))
plot(f_ax,(rad2deg(unwrap(angle(V1dft)))))
% plot(f_ax,(rad2deg(angle(X1dft))))
xlim([0, 1])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Phase Unwrapped (deg)')
ylabel('Phase (deg)')