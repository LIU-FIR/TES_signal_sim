function [y] = pfb_fft(pfb_fir)
%返回pfb_fir的每个通带的平均功率值。
% pfb_fir: PFB FIR Filter，M*N
%    
y = conj(fft(pfb_fir));

end

