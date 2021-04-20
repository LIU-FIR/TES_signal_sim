%%
% Polyphase Filter Bank: Oversampled PFB的Matlab实现
% PFB后加FFT以及PFB后加数字下变频滤波取出所需信号
%%
clear; clc;

M = 4096;% [0,2*pi]划分的子带数目，PFB的分支数，
R = 4;% 每一分支上滤波器的长度(tap)
W = 20;%

L = M*R*W;
N = M*R-1;


%%
% 设计滤波器
Fp = 1/M; % 最大抽取率对应的归一化频率，滤波器带宽
Fsh = 1.6/M; % 过采样率对应的归一化频率
Bf = [-1/M/2 1/M/2]; 
Fc = 1/M; % 截止频率对应的归一化频率
% f = [0 Fp Fst 1];
% a = [1 1 0 0];
% h = firpm(N,f,a);
h = fir1(N,Fc,'low',kaiser(N+1,5));
%%
het1=exp(1i*2*pi*(0:N)*Fsh);
het2=exp(1i*2*pi*(0:N)*(-Fsh));
h1 = h.*het1;
h2 = h.*het2;
r=1;
figure(r);r=r+1;
plot((-0.5:1/16384:.5-1/16384),fftshift(20*log10(abs(fft(h,16384)))),'b')
hold on
plot((-0.5:1/16384:.5-1/16384),fftshift(20*log10(abs(fft(h1,16384)))),'r')
% hold on
% plot((-0.5:1/16384:.5-1/16384),fftshift(20*log10(abs(fft(h2,16384)))),'k')
% xlim([-0.005 0.005])
xlim(20*Bf)
grid on
%%
fs = 4e9;
K = 5; % 设定选定子带的序号，K=0,1,2...
P = 64;
Df = fs/M;
B = 2*pi/M;
% A = 8;
% A = 2;
A = 4;
fA = K*Df+A*Df/P;
wA = 2*pi*fA/fs;
% w0 = K*B; % 第K个子带的中心频率
% % w1 = pi/64+3*B;
% w1 = w0+5*B/64;  % 在第K个子带内设定两个偏移中心频率
% w2 = w0-B/64;
% x = exp(1i*w0*(0:L-1))+1.5*exp(1i*w1*(0:L-1))+4*exp(1i*w2*(0:L-1)); % 第K个子带内模拟信号
x = exp(1i*wA*(0:L-1));


% p=1;q=1;
% y_fir = pfb_fir(x,h,M,R,8,5); % oversampled_PFB，过采样率为P/Q
y_fir = pfb_fir(x,h,M,R,4,3); % oversampled_PFB，过采样率为P/Q
%%
tp= 10*R;
y = pfb_fft(y_fir);
figure(1)
stem((0:M-1),abs(y(:,tp)),'s--', 'MarkerSize',5,'MarkerFaceColor','k')
xlim([0 10])
xlabel('channel')
ylabel('channle power: arbitary unit')
grid on

%%
% 选择第K个通道子带，做FFT谱分析

df = B/P/pi; % 数字频率分辨率
chan = K+1;
f = (0:P-1)*df; % 数字频率轴
fshift = (-P/2:P/2-1)*df;% 以0为中心的数字频率轴
y_sp = fft(y(chan,1:P));
mag_y_sp = abs(y_sp);
figure(2)
% stem(fshift,fftshift(mag_y_sp),'^--', 'MarkerSize',3,'MarkerFaceColor','b')
stem(f,(mag_y_sp),'^--', 'MarkerSize',3,'MarkerFaceColor','b')
xlabel('normalized frequency: (\times \pi /rad/sample)')
ylabel('FFT magnitude: arbitary unit')
grid on
% hold on
% plot(fshift,fftshift(mag_y_sp),'*')
%%
% DDC数字下变频
t = 0:size(y_fir,2)-1;
y_ddc = y_fir(K+1,:).*exp(-1i*w1*t*M); %下变频的数字频率由于子带每隔M*Ts出一个值，因此相位变化为：w*M*Ts
y_ddc_sp = fft(y_ddc(1:P));
mag_y_ddc_sp = abs(y_ddc_sp);
figure(3)
% stem(fshift,fftshift(mag_y_ddc_sp),'^--', 'MarkerSize',3,'MarkerFaceColor','r')
stem(f,(mag_y_ddc_sp),'^--', 'MarkerSize',3,'MarkerFaceColor','r')
grid on
xlabel('normalized frequency: (\times \pi /rad/sample)')
ylabel('FFT magnitude: arbitary unit')
%%
% 下变频后的低通滤波器
Fc1 = 1/(2*P); % 截止频率=pi/P
Hf1 = fdesign.lowpass('N,Fc',P,Fc1);
% Hf = fdesign.lowpass('Fp,Fst,Ap,Ast');
Hd3 = design(Hf1,'window','window',{@chebwin,50}, ...
            'systemobject',true);
hfvt1 = fvtool(Hd3,'Color','White');
h_ddc = Hd3.Numerator;
%%
% 经过低通滤波器后的信号
y_ddc_filt = conv(h_ddc,y_ddc);
y_ddc_filt_sp = fft(y_ddc_filt(1:P));
mag_y_ddc_filt_sp = abs(y_ddc_filt_sp);
figure(4)
stem(fshift,fftshift(mag_y_ddc_filt_sp),'b.')
xlabel('normalized frequency: (\times \pi /rad/sample)')
ylabel('FFT magnitude: arbitary unit')
grid on