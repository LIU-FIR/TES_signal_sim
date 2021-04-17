%%
% Polyphase Filter Bank: PFB的Matlab实现
%%
clear; clc;

M = 4096;
% M=8;R = 20;
R = 4;
P = 40;

L = M*R*P;
N = M*R-1;
%%
% 设计滤波器
Fp = 1/M;
Fst = 1.5/M;
Bf = 2/M;
Fc = 1.2/M;
% f = [0 Fp Fst 1];
% a = [1 1 0 0];
% h = firpm(N,f,a);
h = fir1(N,Fc,'low',chebwin(N+1,50));
%%
het1=exp(1i*2*pi*(0:N)*Fst);
het2=exp(1i*2*pi*(0:N)*(-Fst));
het3=exp(1i*2*pi*(0:N)*(2*Fst));
h1 = h.*het1;
h2 = h.*het2;
h3 = h.*het3;
r=1;
figure(r);r=r+1;
plot((-0.5:1/(M*R):.5-1/(M*R)),fftshift(20*log10(abs(fft(h,M*R)))),'b')
hold on
plot((-0.5:1/(M*R):.5-1/(M*R)),fftshift(20*log10(abs(fft(h1,(M*R))))),'r')
% hold on
% plot((-0.5:1/2048:.5-1/2048),fftshift(20*log10(abs(fft(h2,2048)))),'k')
% hold on
% plot((-0.5:1/2048:.5-1/2048),fftshift(20*log10(abs(fft(h3,2048)))),'g')
% grid on
xlim([-10*Bf 10*Bf])
% Ap = 1;
% Astop = 40;
% Hf = fdesign.lowpass('N,Fc',N,Fc);
% % Hf = fdesign.lowpass('Fp,Fst,Ap,Ast');
% Hd2 = design(Hf,'window','window',{@chebwin,50}, ...
%             'systemobject',true);
% % hfvt = fvtool(Hd2,'Color','White');
% h = Hd2.Numerator;




%%

B = 2*pi/M;
w0 = 6*B;
w1 = w0+B/64;
w2 = w0-B/64;
% x = exp(1i*w0*(0:L-1))+1.5*exp(1i*w1*(0:L-1))+4*exp(1i*w2*(0:L-1));
x = exp(1i*w0*(0:L-1));

y_fir = pfb_fir(x,h,M,R,1,1);
% y_fir = pfb_fir_maxd(x,h,M,R);
y_pfb = pfb_fft(y_fir);
% y_os2v1_fir = os2v1_pfb_fir(x,h,M,R);
% y_os2v1_pfb = pfb_fft(y_os2v1_fir);
disp('end here.')
%%
tp= 2*R;
figure(r);r=r+1;
% stem(abs(y_pfb(:,tp)),'r*')
% hold on
% [m,n]=size(y_os2v1_fir);
stem(0:M-1,abs(y_pfb(:,tp)),'b')
xlim([0 10])
% error = abs(y_os2v1_pfb-y_pfb);
% plot(reshape(error,[m*n,1]),'k-');
%%
%PFB的频率响应：在全带宽[0,2*pi]生成多个复正弦信号，以扫频方式画出每个正弦信号的频率响应。
fspan = linspace(0,M*B,8*101);
t=(0:M*R*P-1);
chan = zeros(M,length(fspan));
chan_os = zeros(M,length(fspan));

for k = 1:length(fspan)
    z = exp(1i*fspan(k)*t);   
    Z = pfb_fft(pfb_fir(z,h,M,R,1,1));
    Z2 = pfb_fft(pfb_fir(z,h,M,R,2,1));
    for p = 1:M
        chan(p,k) = abs(Z(p,tp));
        chan_os(p,k) = abs(Z2(p,tp));
    end
end
%%
figure(r);r=r+1;
for kk = 1:2
    plot(fspan/pi, db(chan(kk,:)))
    hold on
end

xlim([0,fspan(end)/pi])
ylim([-150 1])
xlabel('Normalized frequency: (\times \pi /rad/sample)')
ylabel('Channel Magnitude: dB')
grid on

figure(r);r=r+1;
for kk = 1:2
    plot(fspan/pi, db(chan_os(kk,:)))
    hold on
end
% plot(fspan/pi, db(chan(2,:)))
xlim([0,fspan(end)/pi])
ylim([-150 1])
xlabel('Normalized frequency: (\times \pi /rad/sample)')
ylabel('Channel Magnitude: dB')
grid on

%%

