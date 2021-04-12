function [fir_out] = pfb_fir(x,h,M,R,P,Q)
%pfb_fir：多相滤波器的fir滤波器实现
% fir_out: 输出滤波器组输出的数据，[M,R*(W-1)+1]
% x: 输入信号
% h: 原型低通滤波器，长度为M*R
% M：多相滤波器的分支数 (后续FFT的点数，[0,2pi]划分的信道数)
% R：每一分支上滤波器的长度
% P,Q: 过采样率为P/Q，如P=4,Q=3，过采样率为4/3

Ls=M*Q/P;
K = floor((length(x))/Ls); % 
x = x(1:K*Ls); % 将信号截取为处理长度

hp = reshape(h,[M R]); % 按矩阵M*R组织的LPF系数
% regp = zeros(M,R);
reg = zeros(M*R,1);

fir_out = zeros(M,K); % 存储每次load数据时，多相滤波器的输出
fp = 1;
state = 0;
nstate = lcm(M,Ls)/Ls;
for cnt = 1:K %fp: frame pointer, 范围从1到R*W
    reg(Ls+1:end)=reg(1:end-Ls);% 多相寄存器“蛇形”载入数据段
    reg(1:Ls)= x(fp+Ls-1:-1:fp)';%从输入数据2中截取要载入的数据段，同时取共轭
    regp = reshape(reg,[M R]);
    fir_out(:,cnt)=sum(hp.*regp,2);
    fp = fp+Ls;
    
    %循环buffer操作,根据状态交换fir_out列的前M/2点和后M/2点，
    fir_out(:,cnt)=circshift(fir_out(:,cnt),(M-Ls)*state);
    state = mod(state+1,nstate);    
end

end

