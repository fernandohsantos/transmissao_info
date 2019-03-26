clc; close all;clear all;
tic
    
%Informação a ser enviada
tamanho_info = 10^0;    
info(:,1) = round(rand(1,tamanho_info));
%info(:,1) = linspace(1,1,100);%zeros(1,100)';
%info(:,1) = [0 1 0 1];% 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]'; 

%Sincronismo no "0"
%Parâmetros das equações - mestre
a1m = [0.22 0.29]; %sincronizado e não sincronizado
b1m = [0.4 0.4];
c1m =[8.5 8.5];

%Intensidade de acoplamento entre o mestre e o escravo
%eta = linspace(0,1,5);
eta = .3;

%Quantidade de cenários com condições iniciais distintas (Um para cada bit)
cenarios = length(info(:,1));

%Passos da Integração
deltat = 0.10;
npontos = 5000;
t = 0:deltat:(npontos-1)*deltat;

%Valor do trasiente (quantidade de pontos que não serão considerados)
transit = 0;

%Quantidade de SNRs
Nsnrs = 2;
SNRdB = linspace(-20,20,Nsnrs);
SNR = 10.^(SNRdB/10);

%Parâmetros das equações do escravo
a2 = 0.22;
b2 = 0.4;
c2 = 8.5;

%Pontos utilizados para o calculo da diferença entre os patemares
%formados por delta fi
pontos_iniciais = 10; 
pontos_finais = pontos_iniciais;

for valsnr = 1:Nsnrs
    
    %Condições iniciais do mestre
    condicoesiniciaismestre(1,:) = 10*rand(3,1)-5;
    xmestre(1,:) = condicoesiniciaismestre(1,:)';
    
    %Condições iniciais do escravo
    condicoesiniciaisescravo(1,:) = 10*rand(3,1)-5;
    xescravo(1,:) = condicoesiniciaisescravo(1,:)';

        for numcenarios = 1:cenarios, 
        numcenarios;
        
            if mod(numcenarios,100000) == 0;
            numcenarios
            end
       
                    %vetor com os parâmetros, para "despoluir" a entrada da função
                    a1 = a1m(info(numcenarios)+1);
                    b1 = b1m(info(numcenarios)+1);
                    c1 = c1m(info(numcenarios)+1);
                    parametros_mestre = [a1;b1;c1];

                            %Cálculo dos argumentos do Runge-Kutta para o mestre
                            for i = 2:npontos

                                    k1 = mestre(parametros_mestre,xmestre(i-1,:))*deltat;
                                    k2 = mestre(parametros_mestre,xmestre(i-1,:)+ 0.5*k1)*deltat;
                                    k3 = mestre(parametros_mestre,xmestre(i-1,:)+ 0.5*k2)*deltat;
                                    k4 = mestre(parametros_mestre,xmestre(i-1,:)+ k3)*deltat;

                                % Passo para atualizar o novo valor do vetor x  
                            xmestre(i,:) = xmestre(i-1,:) +(k1+2*k2+2*k3+k4)/6;
                            end

        % Armazena os pontos de cada bit em um unico sinal
        %sinal((numcenarios-1)*npontos+1:npontos*numcenarios,1) = xmestre(:,1);
        %Atualiza as condições iniciais do mestre e escravo,
        %utilizando os últimos pontos do mestre e escravo anteriores
        xmestre(1,:) = xmestre(npontos,:);  
                    
                    
        %%  Segunda Parte - Transmite a informação para o escravo;          

        Potsinal = mean(xmestre(:,1).^2);
        Potruido = Potsinal./SNR;        
        
                    %Cria o ruido a ser transmitido junto ao sinal
                    %ruido = sqrt(Potruido(valsnr))*randn(1,npontos);
                    ruido = 0*randn(1,npontos);
                    transmissao = xmestre(:,1) + ruido';  
                    
                    %for numeta = 1:length(eta),
                    %vetor com os parâmetros, para "despoluir" a entrada da função
                    parametros_escravo = [a2;b2;c2;eta];

                            %Cálculo dos argumentos do Runge-Kutta para o escravo
                            for i = 2:npontos

                                    k1 = escravo(parametros_escravo,transmissao(i-1),xescravo(i-1,:))*deltat;
                                    k2 = escravo(parametros_escravo,transmissao(i-1),xescravo(i-1,:)+ 0.5*k1)*deltat;
                                    k3 = escravo(parametros_escravo,transmissao(i-1),xescravo(i-1,:)+ 0.5*k2)*deltat;
                                    k4 = escravo(parametros_escravo,transmissao(i-1),xescravo(i-1,:)+ k3)*deltat;

                            % Passo para atualizar o novo valor do vetor x  
                            xescravo(i,:) = xescravo(i-1,:) +(k1+2*k2+2*k3+k4)/6;
                            end
          %quando não há ruido, ou a potência do mesmo é muito baixa é possível utilizar 
          xescravo(1,:) = xescravo(npontos,:);
          %Caso a potencia do ruido seja alta o sinal não converge para um atrator 
          %xescravo(1,:)= 10*rand(3,1)-5;
                    %end
        %%  Terceira Parte - Calcula as fases e suas diferenças;
        %--------------------------------------------------------------------------------        
        %                        Trasformadas de Hilbert 
        %--------------------------------------------------------------------------------  
        
         xm = hilbert(xmestre(:,1));
         fm = unwrap(angle(xm));
         
         xs = hilbert(xescravo(:,1));
         fs = unwrap(angle(xs));
         
         % Diferença das fases 
         delta = (fm - fs)/pi();
      
 
              delta_inicial = mean(delta(1:pontos_iniciais));
              delta_final = mean(delta(npontos-pontos_finais:npontos));
        
              %diferenca_delta(numcenarios,k) = delta_inicial - delta_final;
               diferenca_delta(numcenarios,valsnr) = abs(delta_inicial - delta_final);
               
%         if mod(numcenarios,10) == 0;
%         numcenarios
%         save('analise_delta.mat','diferenca_delta','info','-v7.3');
%         end
        %sinal((numcenarios-1)*npontos+1:npontos*numcenarios,1) = xmestre(:,1);
        %sinal((numcenarios-1)*npontos+1:npontos*numcenarios,2) = xescravo(:,1);
        end
        
%%
%------------------------------------------------------------------------
%               "Escolhe" o limiar e Cálculo o bit error
%------------------------------------------------------------------------

j=1;
k=1;
    for i=1:length(info)
        %if diferenca_delta_info(i,2) == 0
        if info(i) == 0
        dif_bit_zero(j,valsnr) = diferenca_delta(i,valsnr);
        j = j+1;
        else
        dif_bit_um(k,valsnr) = diferenca_delta(i,valsnr);
        k = k+1;    
        end  
    end
    
%length(dif_bit_zero)+length(dif_bit_um);

limiar(valsnr) = ((max(dif_bit_zero(:,valsnr)))+min(dif_bit_um(:,valsnr)))/2;
%limiar(valsnr) = min(dif_bit_um(:,valsnr));

%limiar = 20;
limiar = linspace(10,10,Nsnrs)
diferenca_delta(:,valsnr) > limiar(valsnr);
info_recebida(:,valsnr) = ans;
info_recebida(:,valsnr) == info(:,1);
acertos(valsnr) = sum(ans);
BER(valsnr)= abs(acertos(valsnr)-length(info))/length(info)
plot(SNRdB,BER)
        %if mod(valsnr,100) == 0;
        %valsnr
        %save('envio_zro_um.mat','BER','SNRdB','limiar','diferenca_delta','info','info_recebida','-v7.3');
        %save('analise_delta.mat','diferenca_delta','info','-v7.3');
        %end
end
toc

%Foi a última que utilizei.


%save('envio_com_ruido5.mat','BER','SNRdB','limiar','diferenca_delta','info','info_recebida','-v7.3')
%valsnr = 1

%'envio_com_ruido está entre -5 a 10
%'envio_com_ruido2 está entre -10 e -5
%'envio_com_ruido3 está em -10 a 10 com 10^4 info
%'envio_com_ruido4 está em -10 a 10 com 10^6 info
%'envio_com_ruido5 está em -20 a 20 com 10^6 info

%limiar_plot = linspace(limiar(:,valsnr),limiar(:,valsnr),length(info));
% pontos_sinal = linspace(1,length(sinal),length(sinal));
% 
%          xmteste = hilbert(sinal(:,1));
%          fmteste = unwrap(angle(xmteste));
%          
%          xsteste = hilbert(sinal(:,2));
%          fsteste = unwrap(angle(xsteste));
%          
%          deltateste = (fmteste - fsteste)/pi();
% 
% figure(1)
% subplot(211)
% plot(pontos_sinal,fmteste,pontos_sinal,fsteste)
% ylabel('\phi_1(t),\phi_2(t)')
% xlabel('n')
% grid on;
% 
% subplot(212)
% plot(pontos_sinal,deltateste)
% ylabel('\frac{\Delta\phi(t)}{\pi}')
% xlabel('n')
% grid on;
%---------------------------------------------------------------------------------

j=1
for j = 1:20
limiar4(j) = 16.2985  
diferenca_delta(:,j) > limiar4(j);
info_recebida(:,j) = ans;
info_recebida(:,j) == info(:,1);
acertos4(j) = sum(ans);
BER4(j)= abs(acertos4(j)-length(info))/length(info)
end

figure(1)
plot(SNRdB,BER,SNRdB,BER2,SNRdB,BER3,SNRdB,BER4)
ylabel('$BER$')
xlabel('$SNRdb$')
grid on;

% 
plot(SNRdB, BER4)





% figure(2)
% plot(pontos_sinal,deltateste)
% ylabel('\frac{\Delta\phi(t)}{\pi}')
% xlabel('n')
% xlim([5000 10000])
% grid on;
% %---------------------------------------------------------------------------------

% figure(3)
% subplot(311)
% plot(info,'o')
% xlabel('bit')
% ylabel('\textbf{B}')
% grid on;
% 
% subplot(312)
% plot(1:length(info),limiar_plot,1:length(info),diferenca_delta(:,valsnr),'o')
% xlabel('bit')
% ylabel('\mu_{AB}')
% grid on;
% 
% subplot(313)
% plot(1:length(info),info_recebida(:,valsnr),'x',1:length(info),info,'o')
% xlabel('bit')
% ylabel('\textbf{\^B}')
% grid on;
%---------------------------------------------------------------------------------

% figure(4)
% plot(1:npontos,xmestre(:,1),1:npontos,xescravo(:,1))
% ylabel('x_1(t), x_2(t)')
% xlabel('npontos')
% grid on;
% 
% figure(5)
% subplot(211)
% plot(info,'o')
% xlabel('bit')
% ylabel('\textbf{B}')
% grid on;
% 
% subplot(212)
% plot(1:length(info),limiar_plot,1:length(info),diferenca_delta(:,valsnr),'o')
% xlabel('bit')
% ylabel('\mu_{AB}')
% grid on;
% 
% %subplot(212)
% figure(2)
% plot(SNRdB,BER,'o')
% grid on;

figure(3)
nbins = 40; 
Npontos = length(limiar);
DeltaX = max(limiar)-min(limiar);
[histo,X]=hist(limiar,nbins);
K = nbins/(DeltaX*Npontos);
histonorm = histo*K 
plot(X,histonorm)
xlabel('L_{AB}')
ylabel('n(L_{AB})')
grid on

% %Criar um limiar
% j = 1;
% k = 1;
% for i=1:length(diferenca_delta)
%         %if diferenca_delta_info(i,2) == 0
%         if info(i) == 0
%         dif_bit_zero(j) = diferenca_delta(i);
%         j = j+1;
%         else
%         dif_bit_um(k) = diferenca_delta(i);
%         k = k+1;    
%         end  
% end
% 
% length(dif_bit_zero)+length(dif_bit_um)
% 
% limiar(cenarios) = (max(dif_bit_zero)+min(dif_bit_um))/2;


% ruido = sqrt(Potruido(k))*randn(1,npontos);
% transmissao = xmestre(:,1)+ ruido';         
%  
% %Caso exista a necessidade de armazenar o sinal
% %Nome_arquivo xmaster xslave
% 
% %%
% Npontos = 1e6;
% Nsnrs = 20;
% sinal = round(rand(1,Npontos));
% Potsinal = mean(sinal.^2);
% 
% SNRdB = linspace(-20,20,Nsnrs);
% SNR = 10.^(SNRdB/10);
% Potruido = Potsinal./SNR;
% 
% for n = 1:Nsnrs,
%     ruido = sqrt(Potruido(n))*randn(1,Npontos);
%     sinalrecebido = sinal+ruido;
%     sinaldecod = (sign(sinalrecebido-1/2)+1)/2;
%     erro = sinal-sinaldecod;
%     numerros = sum(abs(erro));
%     BER(n) = numerros/Npontos;
% end
% 
% hold on;
% plot(SNRdB,BER);
% hold off;