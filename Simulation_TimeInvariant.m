clear
clc
close all
warning off

%% Time-invariant simulation for gOPDC analysis represented in ref [1] (Fig. 3)
%%% Written by: Amir Omidvarnia, 2013

%%% Ref [1]: A.  Omidvarnia,  G.  Azemi,  B.  Boashash,  J.  O.  Toole,  P.  Colditz,  and  S.  Vanhatalo, 
%%% “Measuring  time-varying  information  flow  in  scalp  EEG  signals:  orthogonalized  partial 
%%% directed coherence,”  IEEE  Transactions on Biomedical Engineering, 2013  [Epub ahead of print]

%% Time-invariant MVAR model 
Fs = 1000; % Sampling frequency
Fmax = Fs/2; % Cut off frequency (Hz), should be smaller than Fs/2
Nf = 40;

y = zeros(3000,5);
L = size(y,1); % Number of samples
CH = size(y,2); % Number of channels
a54 = zeros(1,L);
a45 = zeros(1,L);
for n = 4 : L
    y(n,1) = 0.95*sqrt(2)*y(n-1,1) - 0.9025*y(n-2,1) + 10*randn;
    y(n,2) = 0.5*y(n-2,1) + 5*randn;
    y(n,3) = -0.4*y(n-3,1) + 1*randn;
    y(n,4) = -0.5*y(n-2,1) + 0.25*sqrt(2)*y(n-1,4) + 0.25*sqrt(2)*y(n-1,5) + 1.5*randn;
    y(n,5) = -0.25*sqrt(2)*y(n-1,4) + 0.25*sqrt(2)*y(n-1,5) + 2*randn;
end

N_source = 6;
V = rand(CH,N_source);
S = 3*sprand(N_source,L,.5);
y = y + (V*S)';

%% Time-invariant MVAR parameter estimation
[w, A, C, sbc, fpe, th] = arfit(y, 1, 200, 'sbc'); % ---> ARFIT toolbox
[tmp,p_opt] = min(sbc); % Optimum order for the MVAR model

inp_model.data = y;
inp_model.order = p_opt;
PriorMdl = varm(size(y,2),inp_model.order);
EstMdl = estimate(PriorMdl,y);
%% Connectivity measures (PDC, gOPDC etc)
Fmax = Fs/2;
Nf = 40;
[rPDC,GPDC,OPDC,PDC,GOPDC,S] = PDC_dDTF_imag(A,C,p_opt,Fs,Fmax,Nf,EstMdl.Covariance);

%% Plot
%%%%
close all;
figure, % ---> gPDC
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
 
        imagesc(100:L,linspace(0,y_max,Nf),img)
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
figure, % ---> gOPDC
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

        imagesc(100:L,linspace(0,y_max,Nf),img)
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

