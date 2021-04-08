function [reg_out] = pfb_fir(x,h,M,R)
%pfb_fir：多相滤波器的fir滤波器实现
% reg_out: 输出滤波器组输出的数据，[M,R*(W-1)+1]
% x: 输入信号
% h: 原型低通滤波器，长度为M*R
% M：多相滤波器的分支数 (后续FFT的点数，[0,2pi]划分的信道数)
% R：每一分支上滤波器的长度
%   
W = floor(length(x)/(M*R));
x = x(1:M*R*W); % 将信号截取为处理长度
% xp = reshape(x,[M R*W]);
hp = reshape(h,[M R]); % 按矩阵M*R组织的LPF系数
regp = zeros(M,R);
% ind_m = R*(W-1)+1; % load最后一块数据(M*R)时，数据的帧索引
reg_out = zeros(M,R*W); % 存储每次load数据时，多相滤波器的输出
fp = 1;
for cnt = 1:R*W %fp: frame pointer, 范围从1到R*W
    regp(:,2:end) = regp(:,1:end-1);% 数据从左向右载入寄存器
%     regp = xp(:,fp:fp+R-1);
    regp(:,1)=x(fp+M-1:-1:fp)';%从底向上排列，同时取共轭
    reg_out(:,cnt)=sum(hp.*regp,2);
    fp = fp+M;
end

end

