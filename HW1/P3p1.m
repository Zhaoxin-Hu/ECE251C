% Problem 3.1
clc
clf
close all
clearvars

N = 255;
f_ax = [0, 0.25, 0.4, 1];
amp = [1, 1, 0, 0];
w = [1,1];
x = fir2(N,f_ax,amp);
% figure()
% stem([0:N],x)
% title('Type II LPF')
% fvtool(x,'OverlayedAnalysis','phase')

%% Part a
nfft = 1024;
Xdft = fft(x, nfft);
f_ax = linspace(0,2,nfft);
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
% plot(f_ax,(rad2deg(unwrap(angle(Xdft)))))
plot(f_ax,(rad2deg(angle(Xdft))))
xlim([0, 1])
xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Phase Unwrapped (deg)')
ylabel('Phase (deg)')
idx_p = find(f_ax<0.25);
idx_s = intersect(find(f_ax>0.4),find(f_ax<1));
dp = max(abs(abs(Xdft(idx_p))-1));
ds = max(abs(abs(Xdft(idx_s))-0));

%% Part b
x1 = downsample(x,2,0);
nfft = 1024;
X1dft = fft(x1,nfft);
f_ax = linspace(0,2,nfft);
figure()
stem([0:length(x1)-1],x1)
figure()
for i = [1,2]
    ax(i) = subplot(2,1,i);
end
subplot(ax(1))
plot(f_ax,mag2db(abs(X1dft)))
xlim([0, 1])
ylim([-100, 0])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(ax(2))
plot(f_ax,rad2deg(unwrap(angle(X1dft))))
% plot(f_ax,rad2deg(angle(X1dft)))
xlim([0, 1])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Phase Unwrapped (deg)')
% ylabel('Phase (deg)')
idx_p1 = find(f_ax<0.5);
idx_s1 = intersect(find(f_ax>0.8),find(f_ax<1));
dp1 = max(abs(abs(X1dft(idx_p1))-0.5));
ds1 = max(abs(abs(X1dft(idx_s1))-0));

%% Part c
x2 = downsample(x,2,1);
nfft = 1024;
X2dft = fft(x2,nfft);
f_ax = linspace(0,2,nfft);
figure()
stem([0:length(x2)-1],x2)
figure()
for i = [1,2]
    ax(i) = subplot(2,1,i);
end
subplot(ax(1))
plot(f_ax,mag2db(abs(X2dft)))
xlim([0, 1])
ylim([-100, 0])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(ax(2))
plot(f_ax,rad2deg(unwrap(angle(X2dft))))
% plot(f_ax,rad2deg(angle(X2dft)))
xlim([0, 1])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Phase Unwrapped (deg)')
% ylabel('Phase (deg)')
idx_p2 = find(f_ax<0.5);
idx_s2 = intersect(find(f_ax>0.8),find(f_ax<1));
dp2 = max(abs(abs(X2dft(idx_p2))-0.5));
ds2 = max(abs(abs(X2dft(idx_s2))-0));

%% Part e
% h2e = 2*h(4:2:end);
% nfft = 1024;
% H2edft = fft(h2e,nfft);
% fnorm = linspace(0,2,nfft);
% figure()
% for i = [1,2]
%     ax(i) = subplot(2,1,i);
% end
% subplot(ax(1))
% plot(fnorm,mag2db(abs(H2edft)))
% xlim([0, 1])
% ylim([-100, 0])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Magnitude (dB)')
% subplot(ax(2))
% plot(fnorm,rad2deg(unwrap(angle(H2edft))))
% xlim([0, 1])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Phase Unwrapped (\circ)')