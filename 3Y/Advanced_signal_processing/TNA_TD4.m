clearvars; close all; clc;
matlab = 1;

Nord    = 24 ;          % Filter order, an even number
if(matlab)
	hlp = firpm(Nord,[0 .2 .25 1],[1 1 0 0]); % Lowpass filter design
else
	hlp = remez(Nord,[0 .2 .25 1],[1 1 0 0])';
end
D       = Nord/2;       % Propagation delay
hhp     = zeros(size(1:length(hlp))); 
hhp(D+1)= 1;            % Setting the delay path
hhp     = hhp - hlp;    % Complementary operation

Hlp     = fft(hlp, 1024);
Hhp     = fft(hhp, 1024);
nu      = (0:numel(Hlp)-1)/numel(Hlp);

figure;
subplot(211);
    stem(0:Nord', [hlp;hhp]');
    xlabel('samples'); legend('h_{lp}[n]','h_{hp}[n]');
    title('Complementary filters');
subplot(212);
    plot(nu, abs([Hlp;Hhp]));
    xlabel('\nu'); legend('|H_{lp}(\nu)|','|H_{hp}(\nu)|');
    
%% Decomposition of H(\nu) in 2 all-pass filters
[b,a]   = ellip(7,.5,50,.5);
[d1,d2] = tf2ca(b,a);
num     = 0.5*(conv(fliplr(d1),d2) + conv(d1,fliplr(d2)));
den     = conv(d1,d2);
[H,w]   = freqz(b,a,512);
[A01,w] = freqz(num,den,512);
[A0,w]  = freqz(conv(fliplr(d1),d2),den,512);
[A1,w]  = freqz(conv(d1,fliplr(d2)),den,512);

figure;
subplot(311);
    plot(w/(2*pi), abs([H,A01,A0,A1]));
    legend('|H(\nu)|','|A_{01}(\nu)|','|A_0(\nu)|','|A_1(\nu)|'); grid;
    title('Decomposition of H(\nu) in 2 all-pass filters');
subplot(312);
    plot(w/(2*pi), 20*log(abs([H,A01,A0,A1])));
    legend('|H_{dB}(\nu)|','|A_{01dB}(\nu)|','|A_{0dB}(\nu)|','|A_{1dB}(\nu)|'); grid;
subplot(313);
    plot(w/(2*pi), angle([H,A01,A0,A1]));
    legend('Arg(H(\nu))','Arg(A_{01}(\nu))','Arg(A_0(\nu))','Arg(A_1(\nu))'); grid;
    xlabel('\nu'); 
    
% Complementary High-Pass filter H(\nu) in 2 all-pass filters
numc    = 0.5*(conv(fliplr(d1),d2) - conv(d1,fliplr(d2)));
[A01c,w]= freqz(numc,den,512);

figure;
subplot(311);
    plot(w/(2*pi), abs([A01,A01c]));
    legend('|A_{01}(\nu)|','|A_{01c}(\nu)|'); grid;
    title('Complementary filters');
subplot(312);
    plot(w/(2*pi), 20*log(abs([A01,A01c])));
    legend('|A_{01}(\nu)|_{dB}','|A_{01c}(\nu)|_{dB}'); grid;
subplot(313);
    plot(w/(2*pi), angle([A01,A01c]));
    legend('Arg(A_{01}(\nu))','Arg(A_{01c}(\nu))'); grid;
    xlabel('\nu'); 

figure;
subplot(211);
    plot(w/(2*pi), [abs(A01)+abs(A01c), abs(A01.*A01)+abs(A01c.*A01c)]);
    legend('|A_{01}(\nu)|+|A_{01c}(\nu)|','|A_{01}(\nu)|^2+|A_{01c}(\nu)|^2'); grid;
    title('Complementary filters');
subplot(212);
    plot(w/(2*pi), 20*log([abs(A01)+abs(A01c), abs(A01.*A01)+abs(A01c.*A01c)]));
    legend('|A_{01}(\nu)|_{dB}+|A_{01c}(\nu)|_{dB}','|A_{01}(\nu)|^2_{dB}+|A_{01c}(\nu)|^2_{dB}'); grid;
    xlabel('\nu'); 

%% double-complementary (delay-magnitude) low/high-pass filters
Nord    = 18;
blp     = firceqrip(Nord,0.4,[0.04 0.04]); 
fvtool(blp);
bhp     = [zeros(1,Nord/2), 1, zeros(1,Nord/2)]-blp;

Hlp     = fft(blp, 1024);
Hhp     = fft(bhp, 1024);
nu      = (0:numel(Hlp)-1)/numel(Hlp);
Zhp     = roots(bhp);

figure;
subplot(221);
    zplane(blp,1);
    hold on;
    plot(real(Zhp),imag(Zhp),'or');
    hold off;
    legend('Z_{lp}','P_{lp}','unit','Z_{hp}');
subplot(222);
    plot(nu, abs([Hlp;Hhp]));
    xlabel('\nu'); legend('|H_{lp}(\nu)|','|H_{hp}(\nu)|');
subplot(223);
    stem(0:Nord', [blp;bhp]');
    xlabel('samples'); legend('h_{lp}[n]','h_{hp}[n]');
subplot(224);
    plot(nu, abs([Hlp+Hhp]));
    xlabel('\nu'); legend('|(H_{lp}+H_{hp})(\nu)|');

