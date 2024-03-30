% Benjamin Cotte - 03/05/2022
% MF208 Aeroacoustique
% DM2: Analyse des mesures de bruit de profil réalisees dans la soufflerie
% anechoıque de l'Ecole Centrale de Lyon

clear all
close all


% parametres
c = 0.12; % corde du profil (m)
U = 50;   % vitesse de l'ecoulement amont (m/s)
pref = 20e-6; % reference pressure (Pa)

%%
%%% PARTIE 2.1: Analyse des mesures en régime statique (angles d’attaque fixes)
% load microphone data in static regime
% micros_statique_AOAx: signal de pression acoustique (Pa) a l'angle d'attaque x (en degres)
% fs: frequence d'enchantillonnage (Hz)
load('pac_N63_U50_AOA2','micros_statique_AOA2','fs');
load('pac_N63_U50_AOA18','micros_statique_AOA18','fs');
load('pac_N63_U50_AOA30','micros_statique_AOA30','fs');
%%%

%% question 1
echantillon_liste=[micros_statique_AOA2,micros_statique_AOA18,micros_statique_AOA30];
for i =1:3
    echantillon=echantillon_liste(:,i);
    Nsamples=length(echantillon);
    t_statique = Nsamples/fs;
    df = 8; % frequency resolution (Hz)
    Nfft = fs/df; % number of points in the FFT calculation
    Lwindow = 6400; % number of points of the window
    overlap = 0.5; % overlap rate (between 0 and 1)
    noverlap = round(Lwindow*overlap); % number of overlap points
    % on veut tout le signal => on prend tous les points du signal
    window = ones(Lwindow,1); % rectangular window
    % par défaut la routine pwelsh utilise une fenetre de Hamming, on peut
    % donc injecter directement la plage que l'on souhaite étudier. 
    [PSD_pwelch,freq_pwelch] = pwelch(echantillon,window,noverlap,Nfft,fs);
    PSD_pwelch_dB = 10*log10(PSD_pwelch/pref^2);

    figure(21)
    freq_pwelch=freq_pwelch.*c./U;
    semilogx(freq_pwelch(2:end),PSD_pwelch_dB(2:end),'LineWidth',2)
    xlabel('fc / U')
    ylabel('S_{pp} (dB/Hz)')
    xlim([2e-1 3e1])
    ylim([-5 75])
    grid on
    lgd=legend('angle incidence 2°','angle incidence 18°','angle incidence 30°');
    title(lgd,"Angle des différents profils")
    title("DSP des signaux de pression acoustique d'une aile NACA63_3418")
    hold on 
    % on a pas de recouvrement ici d'ou le 0
end 

%% question 2

for i =1:3
    echantillon=echantillon_liste(:,i);
    Nsamples=length(echantillon);
    t_statique = Nsamples/fs;
    df = 8; % frequency resolution (Hz)
    Nfft = fs/df; % number of points in the FFT calculation
    Lwindow = 6400; % number of points of the window
    overlap = 0.5; % overlap rate (between 0 and 1)
    noverlap = round(Lwindow*overlap); % number of overlap points
    % on veut tout le signal => on prend tous les points du signal
    window = ones(Lwindow,1); % rectangular window
    % par défaut la routine pwelsh utilise une fenetre de Hamming, on peut
    % donc injecter directement la plage que l'on souhaite étudier. 
    [PSD_pwelch,freq_pwelch] = pwelch(echantillon,window,noverlap,Nfft,fs);
    PSD_pwelch_dB = 10*log10(PSD_pwelch/pref^2);
    
    % intégration via la méthode des rectangles
    integrale_rect=0;
    dxt=freq_pwelch(2)-freq_pwelch(1);
    f1 = 70; %Hz
    f2 = 1000; %Hz
    %indices valables pour toutes les mesures statiques
    i1 = find(freq_pwelch>=f1,1);
    i2 = find(freq_pwelch>=f2,1);
    % on approxime via la méthode des rectangles
    for j=i1:i2 % on lit les indices dans la variable freq_pwelsh
        integrale_rect=dxt*PSD_pwelch(j); % méthode des rectangles à gauche 
    end
    % resultat
    integrale_rect=10*log10(integrale_rect/(pref*pref)) % on calcule finalement ce qu'on souhaite

    %intégration avec la méthode des trapèzes
    Plage_frequence_integrale = freq_pwelch(i1:i2);
    integrale_trapz = trapz(Plage_frequence_integrale,PSD_pwelch(i1:i2));
    % resultat
    integrale_trapz=10*log10(integrale_trapz/(pref^2)) % on calcule finalement ce qu'on souhaite en dB 
end

%%
%%% PARTIE 2.1: Analyse des mesures en régime statique (angles d’attaque fixes)
% load microphone data in dynamic regime
% micro_dynamique: signal de pression acoustique (Pa)
% k: frequence reduite (voir definition dans l'article)
% alpha0: angle d'attaque moyen (degres)
% alphaM: angle d'attaque maximal au-dessus ou en-dessous de alpha0 (degres)
% fs: frequence d'enchantillonnage (Hz)

%% question 3 tracé de alpha et micro

load('pac_N63_U50_dynamique','micro_dynamique','AOA_dynamique','k','alpha0','alphaM','fs');
t_dynamique = (0:length(micro_dynamique)-1)/fs;
figure(23)
f0=1.33;
plot(t_dynamique.*f0,micro_dynamique,'k.','LineWidth',1)
set(gca,'FontSize',15)
xlabel('f_0 t')
ylabel('p (Pa)')
grid on
title("Signal acoustique dynamique enregistré")
figure(232)
f0=1.33;
plot(t_dynamique.*f0,AOA_dynamique,'k.','LineWidth',1)
set(gca,'FontSize',15)
xlabel('f_0 t ')
ylabel('alpha (degré) ')
grid on
title("Oscillations de l'angle d'attaque de l'aile NACA63_3418")
% on écoute le signal acoustique
sound(micro_dynamique,fs)

%% question 4 tracé de spectrogramme

% Use the spectrogram Matlab routine without overlap
noverlap = 0;
Lwindow=8192;
Nfft=Lwindow;
[~,freq4_specgram,temps_spectro4,DSP_specgram4] = spectrogram(micro_dynamique,Lwindow,noverlap,Nfft,fs);
DSP_specgram4_dB = 10.*log10(DSP_specgram4./(pref^2));
df_specgram4 = fs/Nfft; % frequency resolution
dt_specgram4 = (1-noverlap)*Lwindow/fs; % temporal resolution
figure(24)
pcolor(f0.*temps_spectro4,freq4_specgram.*c./U,DSP_specgram4_dB)
shading flat
colorbar
colormap jet
[mValues , ~] = max(DSP_specgram4_dB);
[mValues , ~] = max(mValues);
caxis([mValues-60 mValues])
ylim([0.1  20]);
set(gca,'FontSize',15)
set(gca, 'YScale', 'log')
xlabel('f_0 t ')
ylabel('f c / U')
title("Spectrogramme pour N=L=8192 sans overlapping")

%% question 5 tracé de spectrogramme

% Use the spectrogram Matlab routine without overlap
noverlap = 0;
Lwindow=3500;
Nfft=8192;
[~,freq5_specgram,temps_spectro5,DSP_specgram5] = spectrogram(micro_dynamique,Lwindow,noverlap,Nfft,fs);
DSP_specgram5_dB = 10*log10(DSP_specgram5./pref^2);
df_specgram5 = fs/Nfft; % frequency resolution
dt_specgram5 = (1-noverlap)*Lwindow/fs; % temporal resolution
figure(25)
pcolor(temps_spectro5.*f0,freq5_specgram.*c./U,DSP_specgram5_dB)
shading flat
colorbar
colormap jet
% on récupère le max du max
[mValues , ~] = max(DSP_specgram5_dB);
[mValues , ~] = max(mValues);
caxis([mValues-60 mValues])
set(gca, 'YScale', 'log')
set(gca,'FontSize',15)
ylim([0.1  20]);
xlabel('f_0 t ')
ylabel('f c / U')
title("Spectrogramme pour N=8192 > L=3500, technique du 0-padding, sans overlapping")

%% question 6 tracé de spectrogramme comme dans l'article mêmes paramètres 

% Use the spectrogram Matlab routine with overlap of 80%
Lwindow=3500;
Nfft=8192;
overlap = 0.8; % overlap rate (between 0 and 1)
noverlap = round(Lwindow*overlap); % number of overlap points
[~,freq6_specgram,temps_spectro6,DSP_specgram6] = spectrogram(micro_dynamique,Lwindow,noverlap,Nfft,fs);
DSP_specgram6_dB = 10*log10(DSP_specgram6/pref^2);
df_specgram6 = fs/Nfft; % frequency resolution
dt_specgram6 = (1-overlap)*Lwindow/fs; % temporal resolution
figure(26)
pcolor(temps_spectro6.*f0,freq6_specgram.*c./U,DSP_specgram6_dB)
shading flat
colorbar
colormap jet
% on récupère le max du max
[mValues , ~] = max(DSP_specgram6_dB);
[mValues , ~] = max(mValues);
caxis([mValues-60 mValues])
set(gca,'FontSize',15)
set(gca, 'YScale', 'log')
xlabel('f_0 t ')
ylabel('f c / U')
ylim([0.1  20]);
title("Spectrogramme pour N=8192 > L=3500, technique du 0-padding, avec 80% d'overlapping")

%% question 7 OASPL

integrale_rect=zeros(1,271); % on mets 271 car c'est la taille de temps_spectro6
dxt=freq6_specgram(2)-freq6_specgram(1); % pas d'intégration des fréquences
% on approxime via la méthode des rectangles
f1 = 70; %Hz
f2 = 1000; %Hz
%indices des bornes de l'intégrale
i1 = find(freq6_specgram>=f1,1);
i2 = find(freq6_specgram>=f2,1);
for i=i1:i2
    integrale_rect=integrale_rect+dxt*DSP_specgram6(i,:); % méthode des rectangles à gauche 
end
OAt=10*log10(integrale_rect/(pref*pref)); % on calcule finalement ce qu'on souhaite

figure(27)
plot(f0*temps_spectro6,OAt,'k','LineWidth',1.3)
xlim([0,1]);
xlabel('f_0 t')
ylabel('OASPL (dB)')
ylim([60 85])
legend('k=0.01')
title("Bruit total OASPL(t)")
figure(272)
plot(f0*temps_spectro6,OAt,'k','LineWidth',1.3)
xlabel('f_0 t')
ylabel('OASPL (dB)')
ylim([60 85])
legend('k=0.01');
title("Bruit total OASPL(t)")

% on essaye de réaliser un moyennage de phase 
dd=round(length(OAt)/5);
moyenne_de_phase=OAt(1,1:dd)+OAt(1,dd+1:dd*2)+OAt(1,dd*2+1:dd*3)+OAt(1,dd*3+1:dd*4)+OAt(1,dd*4+1:dd*5);
moyenne_de_phase=moyenne_de_phase./5;

figure(273)
plot(f0*temps_spectro6(1:dd),moyenne_de_phase,'k','LineWidth',1.3)
xlabel('f_0 t')
ylabel('OASPL (dB)')
ylim([60 85])
legend('k=0.01');
title("Bruit total OASPL(t) moyenné en phase")
