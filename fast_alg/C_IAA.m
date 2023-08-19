function [P, s_hat]=fun_caiaa_power(y_noisy,Mt,Mr, Iter_no)

Kr=10*Mr;
% rad=0:2*pi/Kr:(2*pi/Kr)*(Kr-1);
rad1=linspace(-pi,pi,Kr+1);
rad1(end)=[];
Y_noisy=reshape(y_noisy,Mr,Mt);
Ar=exp(1j*(0:Mr-1).'*rad1);
% Ar=[Ar eye(Mt,Mt)];

%% First step
% tic
s_hat     = Ar'*Y_noisy./Mr; % DAS estimate
P             = sum(s_hat.*conj(s_hat),2)/Mt;
for jj        = 1:Iter_no
    R         = Ar.*P.'*Ar';
    R_invA    = R\Ar;
    part1     = R_invA'*Y_noisy;
    part2     = sum(conj(R_invA).*Ar, 1).';
    s_hat     = part1./part2;
    P1        = sum(s_hat.*conj(s_hat),2)/Mt;
end
% toc

%% Second step
Kt=10*Mt*Mr;
rad=linspace(-pi,pi,Kt+1);
rad(end)=[];
Ar=exp(1j*(0:Mr-1).'*rad);

% tic
R_invAr    = R\Ar;
part1     = R_invAr'*Y_noisy;
part2     = sum(conj(R_invAr).*Ar, 1).';
% toc;
Y_At     = part1./part2;
P=mean(abs(Y_At).^2,2);
L=Kt/Mr;  
offset=floor(L*0.3);
% toc

%% µÚÒ»¶ÎÆµÆ×

index=[Kt-offset+1:Kt 1:L-offset];
rad_tmp=rad(index);
At=exp(1j*Mr*(0:Mt-1).'*rad_tmp);
P_tmp=P(index);
Y_At_tmp=Y_At(index,:);
% tic;
[P_tmp, s_hat]=fun_ciaa_part_power(Y_At_tmp, At,10,P_tmp);
P(1:L-2*offset)=P_tmp(offset+1:L-offset);
flag=L-2*offset;
% toc
%% ÖÐ¼ä¶ÎÆµÆ×
while flag-offset+L <= Kt
    index=flag-offset+1:flag-offset+L;
    rad_tmp=rad(index);
    At=exp(1j*Mr*(0:Mt-1).'*rad_tmp);
    P_tmp=P(index);
    Y_At_tmp=Y_At(index,:);
    [P_tmp, s_hat]=fun_ciaa_part_power(Y_At_tmp, At,10,P_tmp);
    P(flag+1:flag+L-2*offset)=P_tmp(offset+1:L-offset);
%     flag+1:flag+L-2*offset
    flag=flag+L-2*offset;
end
%% ×îÄ©¶ÎÆµÆ×
index=[flag-offset+1:Kt 1:L-Kt+(flag-offset)];
rad_tmp=rad(index);
At=exp(1j*Mr*(0:Mt-1).'*rad_tmp);
P_tmp=P(index);
Y_At_tmp=Y_At(index,:);
[P_tmp, s_hat]=fun_ciaa_part_power(Y_At_tmp, At,10,P_tmp);
P(flag+1:Kt)=P_tmp(offset+1:offset+(Kt-flag));

end