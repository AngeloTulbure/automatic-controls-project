% Progetto a1
% 
% Controlli Automatici T
% Tulbure Angelo Maximilian (0000923220)
% 07/2022
%

clear all, close all,  clc;

omega_plot_min = 1e-3;
omega_plot_max = 1e6;

% variabili utili per mostrare solo i grafici che ci interessano
show_TF = 1;
show_G_estesa = 1;
show_LL = 1;
show_1t = 1;
show_sens = 1;

% DEFINIZIONE PARAMETRI DEL PROGETTO
rs = 1.2;       % tasso di riproduzione cellule suscettibili
rr = 1.1;       % tasso di riproduzione cellule resistenti
k = 100;        % numero massimo di cellule 
gamma = 0.1;    % cambio da resistenti a suscettibili
beta = 0.1;     % cambio da suscettibili a resistenti
alpha = 0.1;    % cambio da suscettibili a resistenti dopo trattamento farmacologico
ms = 0.8;       % mortalità cellule suscettibili
mr = 0.1;       % mortalità cellule resistenti
nse = 50;       % numero di cellule suscettibili
nre = 50;       % numero di cellule resistenti

syms ns nr cf

%% EQUAZIONI DI STATO E DI USCITA
ns_punto = rs*( 1 - (ns + nr)/k )*ns - ms*cf*ns - beta*ns + gamma*nr - alpha*cf*ns;
nr_punto = rr*(1 - (ns + nr)/k )*nr - mr*cf*nr + beta*ns - gamma*nr + alpha*cf*ns;
y = ns; 

%% COPPIA DI EQUILIBRIO
%dobbiamo calcolarci solo la cfe in quanto la traccia ci fornisce già la
%nse e la nre, entrambe pari a 50
cfe = (rs*(1-(nse+nre)/k)*nse)/(ms*nse+alpha*nse) - (beta*nse)/(ms*nse+alpha*nse) + (gamma*nre)/(ms*nse+alpha*nse);   % usiamo una delle due equazioni a piacere

a11 = rs - ((2*rs*nse)/k) - ((rs*nre)/k) - ms*cfe - beta - (alpha*cfe);
a12 = (-rs*nse)/k + gamma;
a21 = (-rr*nre)/k + beta + (alpha*cfe); 
a22 = rr - ((rr*nse)/k) - ((2*rr*nre)/k) - (mr*cfe) - gamma;

b1 = -ms*nse - alpha*nse;
b2 = -mr*nre + alpha*nse;

% matrici del sistema linearizzato attorno all'equilibrio
A = [a11, a12; a21, a22];
B = [b1; b2];
C = [1,0];
D = 0;

%% FUNZIONE DI TRASFERIMENTO
s = tf('s');
sys = ss(A,B,C,D);  %state system, crea il sistema linearizzato
GG = tf(sys);       %transfer function, fa la funzione di trasferimenot (vedere calcoli a mano)
GG_numeratore = GG.Numerator;
%GG = zpk(GG)

%% DIAGRAMMA DI BODE 
if show_TF
    figure(1);
    bode(GG,{omega_plot_min,omega_plot_max});   
    title('Diagramma di Bode');
    grid on, zoom on;
    
    figure(2);
    margin(GG);
    title("Diagramma di Bode con margine di fase");
    grid on, zoom on;
end

%% SPECIFICHE
% REGOLATORE STATICO
% ERRORE A REGIME NULLO
mus = 1.1;  %guadagno del regolatore statico. è un parametro libero. lo inizializziamo a 1.1
R_s = -mus/s;   %abbiamo aggiunto un polo nell'origine perchè non presente
% Guadagno negativo in quanto la GG è negativa

G_e = R_s*GG;   

% SOVRAELONGAZIONE PERCENTUALE
Spercentuale = 0.1;        %Sovraelongazione del 10%
xi = sqrt(log(Spercentuale)^2/(pi^2+log(Spercentuale)^2));     % 0.5912 
Mf = xi*100;    % 59.1155 risulta piu stringente del Mf minima pari a 30 gradi


% TEMPO DI ASSESTAMENTO
Ta1 = 0.4;              %tempo di assestamento all'1% è di 0.4s
omegac_min = 460/(Mf*Ta1);     % 19.4534

% DISTURBO D'USCITA
%disturbo sull uscita d(t) banda range [0,0.075], attenuazione di 55 dB
Ad = 55;                 %abbattimento di 55 dB
omega_d_min = 0.0001;    % lower bound per il plot, idealmente sarebbe 0
omega_d_MAX = 0.075;
%l'intervallo di omega_d sarebbe da 0 a 0.075, lo poniamo in questo modo per facilità di rappresentazione

% DISTURBO DI MISURA
%rumore di misura n(t) banda range[1.5*10^4, 10^6], attenuazione di 80 dB
An = 80;           %abbattimento di 80 dB
omega_n_min = 1.5*10^4;
omega_n_MAX = 10^6;

omegac_MAX = omega_n_min;

%% MAPPATURA SPECIFICHE
if show_G_estesa
    figure(3);
    % DIAGRAMMI DI BODE DI Ge CON SPECIFICHE
    %
    % Dall'HELP di Matlab:
    % [Mag,phase,w]=bode(GG,{omega_plot_min,omega_plot_max});
    [Mag, phase, omega] = bode(G_e, {omega_plot_min, omega_plot_max}, 'k');
    
    % Disturbo d'uscita
    patch([omega_d_min, omega_d_MAX, omega_d_MAX, omega_d_min], [Ad, Ad,-200,-200], 'y', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 

    % Rumore di misura
    patch([omega_n_min, omega_n_MAX, omega_n_MAX, omega_n_min], [-An, -An, 200, 200], 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
    
    % Pulsazione critica
    patch([omega_plot_min, omegac_min, omegac_min, omega_plot_min], [0,0, -300, -300], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
    legend("A_d", "A_n", "\omega_{c,min}");
    grid on, zoom on, hold on;
    
    margin(Mag, phase, omega);
    
    % Margine di fase
    patch([omegac_min, omegac_MAX, omegac_MAX, omegac_min], [-180 + Mf, -180 + Mf, -180, -180], 'g', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
    legend("G_e(j\omega)", "M_f");
    
    grid on, zoom on, hold on;
    
    title("Mappatura specifiche");
end

%% REGOLATORE DINAMICO: RETE ANTICIPATRICE
% SCENARIO B (il vincolo su Mf non è mai rispettato tra omega_c_min e omega_c_max).
% Sfruttiamo una rete anticipatrice
if show_LL
    figure(4);
  
    % Rete anticipatrice
    % omega_c_star e Mf_star dalle specifiche
    Mf_star = Mf;     % Mf_star = 59.1155+0.1
    omega_c_star = omegac_min+30;        %omegac_min = 19.4534 
    [mag_omega_c_star, arg_omega_c_star, omega_c_star] = bode(G_e, omega_c_star);
    

    % M_star = 1/mag_omega_c_star;
    mag_omega_c_star_dB = 20*log10(mag_omega_c_star);
    M_star = 10^(-mag_omega_c_star_dB/20);        %29.667
    phi_star = Mf_star - 180 - arg_omega_c_star;


    % Formule di inversione per calcolare alpha e tau
    phi_star_deg = phi_star*pi/180;
    tau = (M_star - cos(phi_star_deg))/omega_c_star/sin(phi_star_deg);
    alpha_tau = (cos(phi_star_deg) - inv(M_star))/omega_c_star/sin(phi_star_deg);
    alpha = alpha_tau / tau;
    
    % Check cos phi* > 1/M
    check_flag = cos(phi_star_deg) - inv(M_star);
        if check_flag < 0
            disp('Errore: alpha negativo');
            return;
        end
    
    N_Rd = (1 + tau*s); 
    D_Rd = (1 + alpha*tau*s);
    R_d = N_Rd/D_Rd;
    
    RR = R_s*R_d;     % Il regolatore "complessivo" è pari al prodotto del reg. statico e di quello dinamico
    numR = RR.Numerator;
    denomR = RR.Denominator;
    
    LL = RR*GG;   %oppure Rd*G_e dato che => G_e=Rs*GG
    LL = zpk(LL);
    
    [Mag,phase, omega] = bode(LL, {omega_plot_min, omega_plot_max}, 'k');
    
    % Specifiche sul disturbo in uscita
    patch([omega_d_min, omega_d_MAX, omega_d_MAX, omega_d_min], [-200,-200, Ad, Ad], 'y', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
    
    % Specifiche sul rumore di misura
    patch([omega_n_min, omega_n_MAX, omega_n_MAX, omega_n_min], [-An, -An, 200, 200], 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
    
    % Specifiche sulla pulsazione critica
    patch([omega_plot_min, omegac_min, omegac_min, omega_plot_min], [0,0, -300, -300], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
    legend("A_d", "A_n", "\omega_{c,min}");
    
    grid on, zoom on, hold on;
    
    margin(Mag, phase, omega);
    
    % Margine di fase
    patch([omegac_min, omegac_MAX, omegac_MAX, omegac_min], [-180 + Mf, -180 + Mf, -180, -180], 'g', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3); 
    legend("G_e(j\omega)", "M_f");
    
    grid on, zoom on, hold on;
end

%% Prestazioni in risposta al gradino
if show_1t
    F = LL/(1+LL);    % LL = RR*GG

    figure(5);
    WW = -2;   % ampiezza del gradino
    T_simulation = 1;   %simulazione grafico per 1 secondo
    [y_step,t_step] = step(WW*F, T_simulation);
    plot(t_step,y_step,'b');
    title("Prestazioni in risposta al gradino");
    grid on, zoom on, hold on;
    
    % vincolo sovraelongazione
    patch([0,T_simulation,T_simulation,0],[WW*(1+Spercentuale),WW*(1+Spercentuale),WW*2,WW*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
    ylim([WW*1.5,0]);
    
    % vincolo tempo di assestamento all'1%
    e_star = 0.01;  %errore a regime pari a 0.01 (1%) come da specifica
    patch([Ta1,T_simulation,T_simulation,Ta1],[WW*(1-e_star),WW*(1-e_star),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
    patch([Ta1,T_simulation,T_simulation,Ta1],[WW*(1+e_star),WW*(1+e_star),WW*1.5,WW*1.5],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);
    
    
    Legend_step = ["F(s)"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
    legend(Legend_step);
    
    
    %% Prestazioni sul distrurbo in uscita
    figure(6);
    omega_d = 0.01875;
    amplitude_d = 0.3;

    % Funzione di sensitività
    S = 1/(1+LL);

    tt = (0:1e-2:1e3)';
    d1 = amplitude_d*sin(omega_d*tt);
    d2 = amplitude_d*sin(omega_d*tt*2);
    d3 = amplitude_d*sin(omega_d*tt*3);
    d4 = amplitude_d*sin(omega_d*tt*4);
    dd = d1+d2+d3+d4;
    y_d = lsim(S,dd,tt);
    
    hold on, grid on, zoom on
    plot(tt,dd,'g');
    plot(tt,y_d,'b');
    title("Comportamento disturbo di uscita");
    grid on
    legend('dd','y_d');
    
    % calculate gain
    pdd = findpeaks(dd);
    pyd = findpeaks(y_d);
    sprintf("Attenuazione disturbo sull'uscita d(t) %f dB",20*log10(max(pdd)/max(pyd)))
    
%% Prestazioni sul disturbo di misura
    figure(7);
    omega_n = 1.5e4;
    amplitude_n = 0.2;

    %intervallo di campionamento
    tt = (0:1e-7:1e-2)';
    n1 = amplitude_n*sin(omega_n*tt);
    n2 = amplitude_n*sin(omega_n*tt*2);
    n3 = amplitude_n*sin(omega_n*tt*3);
    n4 = amplitude_n*sin(omega_n*tt*4);
    nn = n1+n2+n3+n4;
    y_n = lsim(-F,nn,tt);
    hold on, grid on, zoom on
    plot(tt,nn,'m');
    plot(tt,y_n,'b');
    title("Comportamento disturbo di misura");
    grid on
    legend('nn','y_n');
    
    % calculate gain
    pnn = findpeaks(nn);
    pyn = findpeaks(y_n);
    sprintf("Attenuazione rumore di misura n(t) %f dB",20*log10(max(pnn)/max(pyn)))
end

%% Tutte le funzioni di sensitività
if show_sens
    F = LL/(1+LL);     % LL = RR*GG
    S = 1/(1+LL);
    Q = RR/(1+LL);

    bode_min = 1e-4;
    bode_max = 1e5;

    figure(8);
    hold on;grid on;
    legend(["F(j\omega)"; "S(j\omega)"; "Q(j\omega)"]);
   
    margin(F,{bode_min,bode_max});
    margin(S,{bode_min,bode_max});
    margin(Q,{bode_min,bode_max});
    legend(["F(j\omega)"; "S(j\omega)"; "Q(j\omega)"]);
    title("Funzioni di sensitività relative a L(j\omega)");
end


%variabili da dare a Simulink per congruenza
x1e = nse;
x2e = nre;
ue = cfe;








