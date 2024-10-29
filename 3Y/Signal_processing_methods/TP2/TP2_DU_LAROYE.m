%TP2 Signal avancé SIA | DU | LAROYE
%%
close all ;
clear;
clc ;
format short g ;

%% I) Implémentation de l'algorithme LMS
% Q2
N = 650;  % Longueur du signal
h = [1 0.3 -0.1 0.2]';  % Réponse impulsionnelle du filtre h

% Signal x[n] (bruit blanc)
x = randn(N, 1);

% Signal d[n] en filtrant x[n] par le filtre h
d = filter(h, 1, x);

% Affichage des premiers échantillons de x et d
figure;
subplot(2, 1, 1);
plot(x);
title('Signal x[n] (Bruit blanc)');
xlabel('Échantillons');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(d);
title('Signal d[n] (Signal filtré)');
xlabel('Échantillons');
ylabel('Amplitude');

%% Q3

% Check LMS algorithm function section at the end of the script

%% Q4
% Parameters
P = 40;    % Ordre du filtre
mu = 0.01; % écart-type

[w, y, e] = algolms(x, d, P, mu);

% Abscisse
n = 1:N;  

figure;
plot(n, d, 'b', 'DisplayName', 'Signal désiré d[n]'); hold on;
plot(n, y, 'r', 'DisplayName', 'Signal estimé y[n]');
plot(n, e, 'g', 'DisplayName', 'Erreur e[n]');
legend show;
xlabel('Temps (n)');
ylabel('Amplitude');
title('Comparaison du signal désiré, du signal estimé et de l erreur');
hold off;

%% Q4.3

[w, y, e, w_hist] = algolms(x, d, P, mu);

% Abscisse
n = 1:N;

figure;
for i = 1:4
    subplot(4,1,i);
    plot(n, w_hist(i, :), 'r', 'DisplayName', ['Coefficient estimé w' num2str(i)]); hold on;
    yline(h(i), 'b--', 'DisplayName', ['Vrai coefficient h' num2str(i)]);
    legend show;
    xlabel('Temps (n)');
    ylabel(['w' num2str(i)]);
    title(['Coefficient w' num2str(i) ' estimé vs h' num2str(i)]);
end

%% Q5 Test de l'ordre P

% Parameters
P_values = [5, 10, 20];  % Valeurs de P (ordre du filtre)
noise_level = 0.01;  % Niveau de bruit ajouté au signal désiré

% Test pour les différentes valeurs de P
for P = P_values
    % Générer le filtre FIR de fréquence 0.5 et d'ordre P
    h = fir1(P-1, 0.5);  % Filtre FIR ordre P, fréquence 0.5
    
    d_clean = filter(h, 1, x);  % Signal désiré sans bruit
    d = d_clean + noise_level * randn(N, 1);  % Ajout du bruit
    
    % Algorithme LMS
    [w, y, e, w_hist] = algolms(x, d, P, mu);
    
    % Plot pour chaque P
    figure;
    subplot(2, 1, 1);
    plot(1:N, d, 'b', 'DisplayName', 'Signal désiré d[n]'); hold on;
    plot(1:N, y, 'r', 'DisplayName', 'Signal estimé y[n]');
    xlabel('Temps (n)');
    ylabel('Amplitude');
    title(['Effet de P = ' num2str(P) ' sur le signal désiré vs estimé']);
    legend show;
    
    subplot(2, 1, 2);
    plot(1:N, e, 'g', 'DisplayName', 'Erreur e[n]');
    xlabel('Temps (n)');
    ylabel('Amplitude');
    title(['Effet de P = ' num2str(P) ' sur l''erreur e[n]']);
    legend show;
end

%% Q5 Test de l'écart type mu

% Script pour tester l'effet de mu sur l'algorithme LMS

% Parameters
P = 10;            % P = 10 fixé
mu_values = [0.01, 0.1, 0.5];  % Valeurs de mu à tester
noise_level = 0.01; % Ajout du bruit

% Générer le filtre FIR de fréquence 0.5 et d'ordre P
h = fir1(P-1, 0.5);

% Générer le signal désiré d en filtrant x avec h et ajouter du bruit
d_clean = filter(h, 1, x);  % Signal désiré sans bruit
d = d_clean + noise_level * randn(N, 1);  % Ajout du bruit

% Boucle sur les différentes valeurs de mu
for mu = mu_values

    [w, y, e, w_hist] = algolms(x, d, P, mu);
    
    % Plot
    figure;
    subplot(2, 1, 1);
    plot(1:N, d, 'b', 'DisplayName', 'Signal désiré d[n]'); hold on;
    plot(1:N, y, 'r', 'DisplayName', 'Signal estimé y[n]');
    xlabel('Temps (n)');
    ylabel('Amplitude');
    title(['Effet de \mu = ' num2str(mu) ' sur le signal désiré vs estimé']);
    legend show;
    
    subplot(2, 1, 2);
    plot(1:N, e, 'g', 'DisplayName', 'Erreur e[n]');
    xlabel('Temps (n)');
    ylabel('Amplitude');
    title(['Effet de \mu = ' num2str(mu) ' sur l''erreur e[n]']);
    legend show;
end


%% II) Application

%% 1) Signal audio avec une voix

% Charger le signal audio dans le fichier Voix1.wav
[x, fs] = audioread('Voix1.wav'); % Lire le signal audio

% Écouter le signal audio
%sound(x, fs);

% Charger Rep.dat
Rep_struct = load('Rep.dat', '-mat');
Rep = Rep_struct.RI;

% Afficher la réponse impulsionnelle
figure;
plot(Rep);
title('Réponse impulsionnelle de la chambre');
xlabel('Échantillons');
ylabel('Amplitude');

%%
% Filtrer le signal audio avec la réponse impulsionnelle (convolution)
y_filtered = conv(x, Rep);

% Écouter le signal filtré (avec écho)
%sound(y_filtered, fs);

% Tracer le signal audio filtré
figure;
subplot(2, 1, 1);
plot(x);
title('Signal original Voix1.wav');
xlabel('Échantillons');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(y_filtered);
title('Signal filtré avec réponse impulsionnelle (avec écho)');
xlabel('Échantillons');
ylabel('Amplitude');

%%
% Ajouter un bruit blanc
noise_level = 0.01; % Niveau de bruit ajustable
y_noisy = y_filtered + noise_level * randn(size(y_filtered));

% Écouter le signal filtré plus bruit
%sound(y_noisy, fs);

% Plot signal filtré avec bruit
figure;
subplot(2, 1, 1);
plot(y_filtered);
title('Signal filtré avec écho (sans bruit)');
xlabel('Échantillons');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(y_noisy);
title('Signal filtré avec écho + bruit blanc');
xlabel('Échantillons');
ylabel('Amplitude');


%%
% Parameters
P = 10;    % Ordre du filtre
mu = 0.2;  % Facteur d'apprentissage

[w, y_lms, e] = algolms(x, y_noisy, P, mu);

% Plot
n = P:length(x);  % On commence à partir de l'indice P pour correspondre à la taille de y_lms et e
figure;
subplot(3, 1, 1);
plot(n, y_noisy(n), 'b', 'DisplayName', 'Signal filtré avec bruit');
title('Signal filtré avec écho et bruit blanc');
xlabel('Temps (n)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(n, y_lms(n), 'r', 'DisplayName', 'Signal estimé après LMS');
title('Signal estimé après LMS (annulation d''écho)');
xlabel('Temps (n)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(n, e(n), 'g', 'DisplayName', 'Erreur de filtrage');
title('Erreur après annulation d''écho');
xlabel('Temps (n)');
ylabel('Amplitude');

% Signal après annulation de l'écho
sound(y_lms, fs);

%% 2) Signal audio avec deux voix

% Charger les audios
[farspeech, fs] = audioread('s1.wav'); % Lire le signal audio
[nearspeech, fs] = audioread('s2.wav'); % Lire le signal audio

% Écouter les signaux audios
%sound(farspeech, fs);
%sound(nearspeech, fs);

% Filtrer le signal audio éloigné avec la réponse impulsionnelle (convolution)
farspeech_filtered = conv(farspeech, Rep, 'same');

% Écouter le signal filtré (avec écho)
%sound(farspeech_filtered, fs);

% Créer le dialogue entre le signal éloigné filtré et le signal proche, et écouter le dialogue.
dialogue = farspeech_filtered + nearspeech;

% Écouter le dialogue combiné
%sound(dialogue, fs);

% Paramètres de l'algorithme LMS
P = 20;    % Ordre du filtre (peut être ajusté)
mu = 0.1; % Facteur d'apprentissage (peut être ajusté)

% Appliquer l'algorithme LMS pour annuler l'écho
[w, dialogue_lms, e] = algolms(nearspeech, farspeech_filtered, P, mu);

% Écouter le signal après l'annulation de l'écho
sound(dialogue_lms, fs);
pause(length(dialogue_lms)/fs); % Pause pour écouter le signal traité


%% LMS Algorithm function

function [w, y, e, w_hist] = algolms(x, d, P, mu)
    % Input
    % x: signal d'entrée
    % d: signal de sortie désiré
    % P: ordre du filtre
    % mu: écart-type
    
    N = length(x);
    w = zeros(P, 1); 
    y = zeros(N, 1); 
    e = zeros(N, 1);
    w_hist = zeros(P, N);
    
    % Récurrence pour chaque échantillon
    for n = P:N
        % Prendre les P échantillons les plus récents de x
        x_n = x(n:-1:n-P+1);  % Vecteur de longueur P (fenêtre glissante)
        
        % sortie estimée y[n] = w^T * x_n
        y(n) = w' * x_n;
        
        % erreur e[n] = d[n] - y[n]
        e(n) = d(n) - y(n);
        
        % Update des coefficients du filtre w
        w = w + mu * e(n) * x_n;

        % coeff de l'algo
        w_hist(:, n) = w;
    end
end