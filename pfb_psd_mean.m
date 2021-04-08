function [y_psd_mean] = pfb_psd_mean(pfb_fir,ini_num)
%返回pfb_fir的每个通带的平均功率值。
% pfb_fir: PFB FIR Filter，M*N
% ini_num: 参与平均的帧数
%   
y_psd = abs(conj(fft(pfb_fir))).^2;
N = size(y_psd,2);
M = size(y_psd,1);
L = floor(N/ini_num);
y_psd_mean = zeros(M,L);
fp = 1;
% 每L帧做一次平均
for k = 1:L
    y_psd_mean(:,k) = mean(y_psd(:,fp:fp+ini_num-1),2);
    fp = fp + ini_num;
end

end

