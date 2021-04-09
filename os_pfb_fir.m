function [fir_out] = os_pfb_fir(x,h,M,R,Q)
%oversampled_pfb_fir：多相滤波器的fir滤波器实现,8/5倍过采样率
% fir_out: 输出滤波器组输出的数据
% x: 输入信号
% h: 原型低通滤波器，长度为M*R,
% M：多相滤波器的分支数 (后续FFT的点数，[0,2pi]划分的信道数)
% R：每一分支上滤波器的长度
%   
Q = 8/5;% 过采样率为Q
Ls =M/Q;
K = floor((length(x))/Ls); % 

x = x(1:Ls*K); % 原序列截取数据长度为M*Q*K

hp = reshape(h,[M R]); % 按矩阵M*R组织的LPF多相滤波器矩阵                                      

fir_out = zeros(M,K); % 存储每次load数据时，多相滤波器的输出
reg = zeros(M*R,1);

fp=1;
% circ_buffer = zeros(M,1);
state = 0;
nstate = lcm(M,Ls)/Ls;
for cnt = 1:K % load了K次,计数器范围从1到K

    reg(Ls+1:end)=reg(1:end-Ls);% 多相寄存器“蛇形”载入数据段
    reg(1:Ls)= x(fp+Ls-1:-1:fp)';%从输入数据2中截取要载入的数据段，同时取共轭
    fp = fp+Ls; % 更新当前地址
    regp = reshape(reg,[M R]);
    fir_out(:,cnt)=sum(hp.*regp,2); % 卷积计算PFB_fir相应列输出,从第1列开始，(列下标从1开始)
    
    %循环buffer操作,根据状态交换fir_out列的前M/2点和后M/2点，
    fir_out(:,cnt)=circshift(fir_out(:,cnt),(M-Ls)*state);
    state = mod(state+1,nstate);

end
 
end

