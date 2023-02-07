EbNo = -10:5:25;   
N = 8;              % Tx antennas
M = 8;              % Rx antennas
ch = min([M, N]);
sim_times = 1000;

for i=1:length(EbNo)
    snr = EbNo(i)
    for j=1:sim_times
        H = sqrt(1/2)*(randn(ch) + j*randn(ch));
        [U,S,V] = svd(H);
        
        for k =1:ch
            lambda(k) = S(k,k);
        end
        %lambda = lambda.^2;
        W = sqrt(1/2)*(randn(ch) + j*randn(ch));
        No = 1;

        P = No/(10^(-snr/10));
        optimal = water_filling(lambda, No, P);
        C_opt(j) = sum(log2(1 + (lambda.*optimal)*(1/No)));
        C_equ(j) = sum(log2(1 + (lambda*(P/length(lambda)))*(1/No)));
        C_sin(j) = log2(1 + (max(lambda)*P)*(1/No));
        
        
    end
    C_opt_snr(i) = mean(C_opt);
    C_equ_snr(i) = mean(C_equ);
    C_sin_snr(i) = mean(C_sin);
end

figure(1)
plot(EbNo, C_opt_snr, 'DisplayName','C_optimal');
hold on
plot(EbNo, C_equ_snr, 'DisplayName','C_equal');
plot(EbNo, C_sin_snr, 'DisplayName','C_sin');
legend
hold off
    
% function P_alloc = water_filling(lambda, No, P)
%     sum = 0;
%     for i =1:length(lambda)
%         sum = sum + (No/lambda(i));
%     end
%     
%     V = (P+sum)/length(lambda);
%     
%     for i = 1:length(lambda)
%         P_alloc(i) = max([0 V-(No/lambda(i))]);
%     end
%     
% end

function P_alloc = water_filling(lambda, No, P)
    N=length(lambda);
    No_Lambda = No./lambda;
    V=min(No_Lambda) + P/N; 
    Pv=sum(max(V-No_Lambda,0)); 

    while abs(P-Pv) > 1e-5
        V = V + (P-Pv)/N;
        Pv = sum(max(V-No_Lambda,0));
    end
    P_alloc = max(V-No_Lambda,0);
end