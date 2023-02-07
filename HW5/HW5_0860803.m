EbNo = -10:5:25;   
N = 4;              % Tx antennas
M = 4;              % Rx antennas
ch = min([M, N]);
sim_times = 1000;

Q = eye(ch);
H = sqrt(1/2)*(randn(ch) + j*randn(ch));

 
for i=2:sim_times
   H(:,:,i)  = sqrt(1/2)*(randn(ch) + j*randn(ch)); %Matrix Realizations
end

for k=1:ch
    for i=1:sim_times
        QH(i,k) = sum((Q*H(:,k,i)).^2);
    end
end


No = 1;
for i=1:length(EbNo)
    snr = EbNo(i)
    P = No/(10^(-snr/10));
    
%     for k=1:ch 
%         for s=1:sim_times
%             H_ = H(:,:,s);
%             W = inv(((H_')*H_) + (P/ch)*(1/No).*eye(ch))*H_';
%             WH(s,k) = sum((W*H(:,k,s)).^2);
%         end
%     end

    for s=1:sim_times
        H_ = H(:,:,s);
        No = 1/ch;
        for k=1:ch
            K = (No.*eye(ch)) + (P/ch)*Terms_sum(H_, k);
            KH(s,k) = (H_(:,k)')*(inv(K)*H_(:,k));
            K = (No.*eye(ch)) + (P/ch)*Terms_sum_SIC(H_, k);
            KH_SIC(s,k) = (H_(:,k)')*(inv(K)*H_(:,k));
        end 
        % R_MMSE_SIC(s) = log2(det(eye(ch) + H_*(Q*(P/ch))*H_'));
    end

    R_ZF = 0;
    R_MMSE = 0;
    R_MMSE_SIC = 0;
    for c=1:ch
        R_ZF = R_ZF + mean(log2(1 + QH(:,k).*(P/ch)));
        R_MMSE = R_MMSE + mean(log2(1 + KH(:,k).*(P/ch)));
        R_MMSE_SIC = R_MMSE_SIC + mean(log2(1 + KH_SIC(:,k).*(P/ch)));
    end
    R_ZF_snr(i) = R_ZF;
    R_MMSE_snr(i) = R_MMSE;
    R_MMSE_SIC_snr(i) = mean(R_MMSE_SIC);
    
end



figure(1)
plot(EbNo, R_ZF_snr, 'DisplayName','ZF');
hold on
plot(EbNo, R_MMSE_snr, 'DisplayName','MMSE');
plot(EbNo, R_MMSE_SIC_snr, 'DisplayName','MMSE SIC');
title('Rate vs SNRdb')
legend;
hold off

function y = Terms_sum(H_, k)
    temp = H_;
    temp(:,k) = [];
    y = 0;
    for s=1:length(temp)-1
        y = y + temp(:,s)*(temp(:,s)');
    end
end

function y = Terms_sum_SIC(H_, k)
    temp = H_;
    y = 0;
    for s=(k+1):length(temp)
        y = y + temp(:,s)*(temp(:,s)');
    end
end