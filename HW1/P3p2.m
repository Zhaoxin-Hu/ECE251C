% Problem 3.2
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
% stem([0:255],h)
% title('Type II LPF Coeff')
% fvtool(h)

%% Part a
% y = 3*downsample(x,3,0);
% nfft = 1024;
% Ydft = fft(y,nfft);
% f_ax = linspace(0,2,nfft);
% figure()
% stem([0:length(y)-1],y)
% figure()
% for i = [1,2]
%     ax(i) = subplot(2,1,i);
% end
% subplot(ax(1))
% plot(f_ax,mag2db(abs(Ydft)))
% plot(f_ax,abs(Ydft))
% xlim([0, 1])
% % ylim([-10, 0])
% ylim([0, 1])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% % ylabel('Magnitude (dB)')
% ylabel('Magnitude')
% subplot(ax(2))
% plot(f_ax,rad2deg(unwrap(angle(Ydft))))
% % plot(f_ax,rad2deg(angle(X1dft)))
% xlim([0, 1])
% xlabel('Normalized frequency (\times \pi rad/sample)')
% ylabel('Phase Unwrapped (deg)')
% % ylabel('Phase (deg)')
% idx_py = find(f_ax<0.75);
% idx_sy = intersect(find(f_ax>0.8),find(f_ax<1));
% dpy = max(abs(abs(Ydft(idx_py))-1));
% dsy = max(abs(abs(Ydft(idx_sy))-0));

%% Part b
z = 4*downsample(x,4,0);
nfft = 1024;
Zdft = fft(z,nfft);
f_ax = linspace(0,2,nfft);
figure()
stem([0:length(z)-1],z)
figure()
for i = [1,2]
    ax(i) = subplot(2,1,i);
end
subplot(ax(1))
plot(f_ax,mag2db(abs(Zdft)))
% plot(f_ax,abs(Zdft))
xlim([0, 1])
ylim([-10, 10])
% ylim([0, 1])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
% ylabel('Magnitude')
subplot(ax(2))
plot(f_ax,rad2deg(unwrap(angle(Zdft))))
% plot(f_ax,rad2deg(angle(X1dft)))
xlim([0, 1])
xlabel('Normalized frequency (\times \pi rad/sample)')
ylabel('Phase Unwrapped (deg)')
% ylabel('Phase (deg)')
% idx_py = find(f_ax<0.75);
% idx_sy = intersect(find(f_ax>0.8),find(f_ax<1));
% dpy = max(abs(abs(Zdft(idx_py))-1));
% dsy = max(abs(abs(Zdft(idx_sy))-0));