clear
clc
close all
warning off

%% Time-varying simulation for gOPDC adn rPDC analysis 

Fs = 200; % Sampling frequency (there is no difference between different Fs values for the simulated data. PDCs are the same)
L = 1000; % Number of samples
y = zeros(L,3);
CH = size(y,2); % Number of channels
c = zeros(1,L);

bb = sinc(linspace(pi/2+pi/4,5*pi,L));
b = (1*(bb-min(bb))/(max(bb)-min(bb)))-.2; % Time-varying MVAR parameter

for n = 3 : L
    if(n<=L/2)
        c(n) = (n/(L/2));
    else
        c(n) = (L-n)/(L/2); % Time-varying MVAR parameter
    end
    y(n,1) = 0.59*y(n-1,1) - 0.2*y(n-2,1) + b(n)*y(n-1,2) + c(n)*y(n-1,3) + randn;
    y(n,2) = 1.58*y(n-1,2) - 0.96*y(n-2,2) + randn;
    y(n,3) = 0.6*y(n-1,3)  - 0.91*y(n-2,3) + randn;

end
p = 2;

%% Add a linear superposition of K independent sources
N_source = 6;
V = rand(CH,N_source);
S = sprand(N_source,L,.5);
y = y + (V*S)';
save('EEG_test.mat','y')
%% DEKF for time-varying MVAR parameter estimation
inp_model.data = y;
inp_model.order = p;
[A,C] = DEKF3(inp_model);        % Estimated time-varying parameters, A = [A1 A2 ... Ar]
PriorMdl = varm(size(y,2),inp_model.order);
EstMdl = estimate(PriorMdl,y);
%% Connectivity measures (PDC, gOPDC etc)
Fmax = Fs/2;
Nf = 40;
[rPDC,GPDC,OPDC,PDC,GOPDC,S] = PDC_dDTF_imag(A,C,p,Fs,Fmax,Nf,EstMdl.Covariance);

%% Plot
%%%%
close all;
%figure, % ---> gPDC
s1 = 0;
clear mask
rPDC = abs(rPDC);
x_max = L;
y_max = Fmax;

for i = 1 : CH
    for j = 1 : CH
        s1 = s1 + 1;
        h = subplot(CH,CH,s1);
        set(h,'FontSize',20,'FontWeight','bold');
        
        rPDC_tmp = (rPDC(i,j,:,:));
        rPDC_tmp=rPDC_tmp./max(max(rPDC_tmp));
        img = squeeze(rPDC_tmp);
        img = img(:,100:end);
        if(i==j)
            img = zeros(size(img));
        end
 rPDC_mean(i,j)=mean(mean(img));
      %  imagesc(100:L,linspace(0,y_max,Nf),img)
        set(h,'YDir','normal')
        caxis([0 1])
        grid on
        
        if(i<CH || j>1)
            set(h,'YTick',zeros(1,0),'XTick',zeros(1,0));
        end

        if(i==CH && j==ceil(CH/2))
            xlabel('Time (sample)','Fontsize',20,'FontWeight','bold')
        end
        if(i==ceil(CH/2) && j==1)
            ylabel('Frequency (Hz)','Fontsize',20,'FontWeight','bold')
        end
        
    end
end
h2 = colorbar;
set(h2, 'Position', [.92 .11 .03 .8150],'FontSize',20)

%%%%
%figure, % ---> gOPDC
s1 = 0;
clear mask
GOPDC = abs(GOPDC);
for i = 1 : CH
    for j = 1 : CH
        s1 = s1 + 1;
        h = subplot(CH,CH,s1);
        set(h,'FontSize',20,'FontWeight','bold');
        
        GOPDC_tmp = (GOPDC(i,j,:,:));
        GOPDC_tmp=GOPDC_tmp./max(max(GOPDC_tmp));
        img = squeeze(GOPDC_tmp);
        img = img(:,100:end);
        if(i==j)
            img = zeros(size(img));
        end
 gOPDC_mean(i,j)=mean(mean(img));
      %  imagesc(100:L,linspace(0,y_max,Nf),img)
        set(h,'YDir','normal')
        caxis([0 1])
        grid on
        
        if(i<CH || j>1)
            set(h,'YTick',zeros(1,0),'XTick',zeros(1,0));
        end

        if(i==CH && j==ceil(CH/2))
            xlabel('Time (sample)','Fontsize',20,'FontWeight','bold')
        end
        if(i==ceil(CH/2) && j==1)
            ylabel('Frequency (Hz)','Fontsize',20,'FontWeight','bold')
        end
        
    end
end
h2 = colorbar;
set(h2, 'Position', [.92 .11 .03 .8150],'FontSize',20)


gOPDC_mean=gOPDC_mean./max(max(gOPDC_mean))
rPDC_mean=rPDC_mean./max(max(rPDC_mean))


save('Connections.mat','rPDC_mean','gOPDC_mean')