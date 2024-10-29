%% TP3 Signal avancé SIA | DU Jack | LAROYE Ambroise
%%
close all ;
clear;
clc ;
format short g ;

%% I) PART 1
% 1) Rappeler les équations de l'algorithme RLS
%% Q2) Implémentation de l'algorithme RLS

% Check at the bottom of the file, the 'rls_filter' function.

%% Q3) Validation de l’algorithme RLS
% Paramètres
L = 4;               % Taille du filtre RLS
lambda = 0.99;       % Facteur d'oubli
delta = 1;           % Paramètre de régularisation
N = 150;            % Nombre d'échantillons

% Définition du filtre inconnu
h = [1; 0.8; -0.3; 0.4];

% Signal d'entrée (bruit blanc)
x = randn(N, 1);

% Signal de référence (filtrage de x avec le filtre h)
d = filter(h, 1, x);

% Appel fonction RLS
[y, e, w, w_evolution] = rls_filter(x, d, L, lambda, delta);

% Plot
figure;
subplot(2,1,1);
plot(e);
title('Erreur de sortie en fonction du temps');
xlabel('Échantillons'); ylabel('Erreur e(n)');

% Tracé des vrais coefficients et des estimations RLS
subplot(2,1,2);
plot(repmat(h, 1, N)', '--', 'LineWidth', 1.5); hold on;
plot(w_evolution', 'LineWidth', 1);
title('Évolution des coefficients du filtre estimés par RLS');
xlabel('Échantillons'); ylabel('Amplitude des coefficients');
legend('h_1 (vrai)', 'h_2 (vrai)', 'h_3 (vrai)', 'h_4 (vrai)', ...
       'w_1 (RLS)', 'w_2 (RLS)', 'w_3 (RLS)', 'w_4 (RLS)');
grid on;


%% Q4) Test de l’algorithme RLS avec un signal simulé

% Paramètres
N = 500;                % Nombre d'échantillons
delta = 1;              % Paramètre de régularisation
noise_std = 0.1;        % Écart-type du bruit ajouté

% Signal d'entrée (bruit blanc)
x = randn(N, 1);

% Filtre de fréquence de coupure fc = 0.5 (FIR de longueur P)
for P = [5, 10]  % Longueur du filtre
    b = fir1(P-1, 0.5);  % Coefficients du filtre FIR

    % Signal de référence d (filtrage de x plus ajout d'un bruit)
    d = filter(b, 1, x) + noise_std * randn(N, 1);

    % Étude de l'effet de lambda
    for lambda = [0.1, 0.5, 0.9]
        % Application de l'algorithme RLS
        [y, e, w, w_evolution] = rls_filter(x, d, P, lambda, delta);

        % Plot
        figure;
        subplot(3,1,1);
        plot(d); hold on; plot(y);
        title(['Signal de référence d et sortie estimée y, P = ' num2str(P) ', \lambda = ' num2str(lambda)]);
        legend('Signal désiré d', 'Sortie estimée y');
        xlabel('Échantillons'); ylabel('Amplitude');

        subplot(3,1,2);
        plot(e);
        title(['Erreur e(n) en fonction du temps, P = ' num2str(P) ', \lambda = ' num2str(lambda)]);
        xlabel('Échantillons'); ylabel('Erreur e(n)');

        subplot(3,1,3);
        plot(w_evolution');
        title(['Évolution des coefficients du filtre estimés, P = ' num2str(P) ', \lambda = ' num2str(lambda)]);
        xlabel('Échantillons'); ylabel('Amplitude des coefficients');
        grid on;
    end
end

%% PARTIE 2

% Charger le signal audio (Voix1.wav)
[x, fs] = audioread('Voix1.wav');

% Écouter le signal audio
%sound(x, fs);

% Longueur du signal
N = length(x);

% Bruit de référence N1 comme un bruit blanc
N1 = randn(N, 1);  % Bruit blanc non corrélé au signal audio

% Filtre FIR de fréquence 0.5 et ordre 32
ordre = 32;
fc = 0.5;
h = fir1(ordre, fc); % Coefficients du filtre FIR

% Bruit filtré N0, corrélé à N1
N0 = filter(h, 1, N1);

% Ajouter le bruit N0 au signal audio
x_bruite = x + N0;

% Plot
figure;
subplot(2, 1, 1);
plot((0:N-1)/fs, x);
title('Signal audio original');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot((0:N-1)/fs, N1);
title('Bruit de référence N1');
xlabel('Temps (s)');
ylabel('Amplitude');


% Plot
figure;
subplot(3, 1, 1);
plot((0:N-1)/fs, x);
title('Signal audio original');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot((0:N-1)/fs, N0);
title('Bruit corrélé N0 (filtré de N1)');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot((0:N-1)/fs, x_bruite);
title('Signal audio avec bruit N0 ajouté');
xlabel('Temps (s)');
ylabel('Amplitude');

%% Ecouter le signal filtré
%sound(x_bruite, fs);

audiowrite('Voix1_bruite.wav', y, fs);
%%

% Ajouter le bruit N0 au signal audio
d = x + N0;  % Signal bruité

% Paramètres
L = 32;          % Taille du filtre RLS (doit correspondre à l'ordre du filtre FIR)
lambda = 0.99;   % Facteur d'oubli
delta = 1;       % Paramètre de régularisation

% Estimer le signal s(n) en utilisant l'algorithme RLS
[y, e, w, w_evolution] = rls_filter(N1, d, L, lambda, delta);  % N1 comme signal de référence


% Plot
figure;

subplot(3, 1, 1);
plot((0:N-1)/fs, x);
title('Signal audio original s(n)');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot((0:N-1)/fs, y);
title('Signal estimé par RLS');
xlabel('Temps (s)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot((0:N-1)/fs, e);
title('Erreur de l''estimation (s(n) - y(n))');
xlabel('Temps (s)');
ylabel('Amplitude de l''erreur');

%% TEST soustraction adaptative de bruit
sound(e, fs);
audiowrite('Voix1_soustract.wav', e, fs);


%% III. PARTIE 3

%% 2) Mise en œuvre de l algorithme NLMS
% Check below the filter_nlms function

%% 3) Test de l'algorithme NLMS avec un signal simulé

% Paramètres
fs = 1000;         % Fréquence d'échantillonnage
N = 650;           % Nombre d'échantillons
P_values = [5, 10, 20];  % Longueurs du filtre
mu_values = [0.01, 0.2, 0.5];

% Boucle sur les différentes longueurs de filtre P
for P = P_values
    % Générer un bruit blanc
    x = randn(N, 1);  % Signal d'entrée (bruit blanc)
    
    % Créer un filtre FIR de réponse impulsionnelle
    h = fir1(P - 1, 0.5);  % Filtre FIR de longueur P et fréquence de coupure 0.5
    
    % Filtrer le bruit blanc pour obtenir le signal désiré
    d_clean = filter(h, 1, x);  % Signal désiré (sans bruit)
    
    % Ajouter un bruit au signal désiré
    d = d_clean + 0.1 * randn(N, 1);
    % Tester avec différents mu
    for mu = mu_values
        % Estimer le signal en utilisant l'algorithme NLMS
        [y, e, w] = filter_nlms(x, d, P, mu);
        
        % Plot
        figure;

        subplot(3, 1, 1);
        plot(d);
        title(['Signal désiré d(n) avec P = ', num2str(P)]);
        xlabel('Echantillons');
        ylabel('Amplitude');

        subplot(3, 1, 2);
        plot(y);
        title(['Signal estimé par NLMS avec \mu = ', num2str(mu)]);
        xlabel('Echantillons');
        ylabel('Amplitude');

        subplot(3, 1, 3);
        plot(e);
        title('Erreur de l''estimation (d(n) - y(n))');
        xlabel('Echantillons');
        ylabel('Amplitude de l''erreur');
    end
end

%% Test on Voix1.wav

% Paramètres
[signal, fs] = audioread('Voix1.wav');   % Charger le signal audio
N = length(signal);                      % Nombre d'échantillons du signal audio
P = 5;                                   % Longueurs du filtre
mu = 0.5;                                % Valeurs de mu


% Créer un filtre FIR de réponse impulsionnelle
h = fir1(P - 1, 0.5);  % Filtre FIR de longueur P et fréquence de coupure 0.5

% Filtrer le signal audio pour obtenir le signal désiré
d_clean = filter(h, 1, signal);  % Signal désiré (sans bruit)

% Ajouter un bruit au signal désiré
d = d_clean + 0.1 * randn(N, 1);  % Signal bruité


% Estimer le signal en utilisant l'algorithme NLMS
[y, e, w] = filter_nlms(signal, d, P, mu);

% Plot
figure;

subplot(3, 1, 1);
plot(d);
title(['Signal désiré d(n) avec P = ', num2str(P)]);
xlabel('Echantillons');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(y);
title(['Signal estimé par NLMS avec \mu = ', num2str(mu)]);
xlabel('Echantillons');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(e);
title('Erreur de l''estimation (d(n) - y(n))');
xlabel('Echantillons');
ylabel('Amplitude de l''erreur');

%% TEST SOUND
sound(y, fs);


%% Fonction RLS

function [y, e, w, w_evolution] = rls_filter(x, d, L, lambda, delta)

    % Input:
    % x      - Signal d'entrée
    % d      - Signal de référence (ou cible)
    % L      - Taille du filtre
    % lambda - Facteur d'oubli (typiquement proche de 1, ex : 0.99)
    % delta  - Paramètre de régularisation pour initialiser la matrice de covariance
    %
    % Output:
    % y           - Sortie du filtre
    % e           - Erreur entre d et y
    % w           - Coefficients finaux du filtre
    % w_evolution - Évolution des coefficients du filtre au fil du temps

    N = length(d);

    % Initialisation
    w = zeros(L, 1);          % Vecteur des poids du filtre
    P = delta * eye(L);       % Matrice de covariance initiale
    y = zeros(N, 1);          % Sortie du filtre
    e = zeros(N, 1);          % Erreur entre d et y
    w_evolution = zeros(L, N); % Évolution des coefficients estimés

    % Adaptation du filtre RLS
    for n = L:N
        % Vecteur d'entrée
        x_n = x(n:-1:n-L+1);
        
        % Sortie du filtre
        y(n) = w' * x_n;
        
        % Erreur
        e(n) = d(n) - y(n);
        
        % Mise à jour RLS
        Pi = P * x_n;                    % Produit de la matrice de covariance par x_n
        K = Pi / (lambda + x_n' * Pi);   % Gain de Kalman
        w = w + K * e(n);                % Mise à jour des poids
        P = (P - K * x_n' * P) / lambda; % Mise à jour de la matrice de covariance
        
        % Valeurs des coefficients pour observer l'évolution
        w_evolution(:, n) = w;
    end
end

%%

%% NLMS algorithm function

function [y, e, w] = filter_nlms(x, d, L, mu)
    % Args:
    % x - signal d'entrée
    % d - signal désiré
    % L - ordre du filtre
    % mu - taux d'apprentissage
    
    N = length(x);       % Nombre d'échantillons
    w = zeros(L, 1);     % Initialisation des coefficients du filtre
    y = zeros(N, 1);     % Sortie du filtre
    e = zeros(N, 1);     % Erreur

    % Boucle principale
    for n = L:N
        % Vecteur d'entrée du filtre de taille L
        x_n = x(n:-1:n-L+1);
        
        % Sortie du filtre
        y(n) = w' * x_n;

        % Calculer l'erreur
        e(n) = d(n) - y(n);                        % Erreur de prédiction
        
        % Mise à jour des coefficients
        norm_x_n = norm(x_n)^2;                    % Norme du vecteur d'entrée
        if norm_x_n > 0                            % Pour éviter la division par zéro
            w = w + (mu / norm_x_n) * e(n) * x_n;  % Mise à jour des coefficients
        end
    end
end
