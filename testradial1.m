%% Temperature dans un fil électrique, flux nulle au centre et flux imposée à l'extremité
% deux couches dans le câble, au centre du cuivre. puis du caoutchou

%% Parametres Matlab
clear all                        
close all                        
clc                              
format long   
tic

%% Parametres de simulation
R              = 1; % rayon de câble
Rc             = 0.4 % rayon de cuivre
N              = 19; % nombre d'elements de taille dx=L/N
t_final        = 0.8; % duree de la simulation
average        = 'arithmetic'; % evaluation de la conductivite aux demi-elements
                               % par moyenne : arithmetic, harmonic ou geometric 
n_frequency    = 10000; % retour ecran toutes les n_frequency iterations temporelles
mur            = '4'; % Mur, 1:simple,  2:tri-couche, 3:conductivite en loi de puissance, 4:cable double couche
phi            = 2;   % voir  
S              = 5;   % valeur de source dépend de amperage à ajuster/calculer
h              = 10;  % important convection pour flux                               
%% Proprietees du mur
switch mur
    case '1' % Mur simple
        k1   = 1; 
        rho1 = 1;   
        c1   = 1;   
    case '2' % Mur tri-couche
        r12 = 0.3*R; % separation couches 1-2
        r23 = 0.5*R; % separation couches 2-3
        % Couche 1
        k1   = 1000; 
        rho1 = 1;    
        c1   = 1;    
        % Couche 2
        k2   = 1000; 
        rho2 = 1;    
        c2   = 1;    
        % Couche 3
        k3   = 1; 
        rho3 = 1;    
        c3   = 1;    
    case '3' % Mur avec conductivite k en loi de puissance k=b*T^a
        rho1 = 1;    
        c1   = 1;    
        b    = 1;
        a    = 3;
    case '4' %Mur double-couche
        r12 = Rc*R;
        % Couche 1
        kc   = 380; 
        rhoc = 8.96;    
        cc   = 385;    
        % Couche 2
        ki   = 0.16; 
        rhoi = 0.93;    
        ci   = 2092;  
end

%% Condition initiale et conditions limites
Tleft  = 25;   % temperature cote gauche TEMPERATURE AU D2PART DANS LE CABLE (TEMPERATURE EXTERIEUR)
Tright = 25; % temperature cote droit TEMPERATURE A LEXTERIEUR
Tmid   = 25; % temperature initiale dans le reste du domaine

%% Espace
dr = R/N;
drm2 = 1/dr/dr;
r = linspace(0,R,N+1);

%% Estimation de alpha_c un coef de diffusion caracteristique du probleme
switch mur
    case '1' % Mur simple, conductivite constante
        alpha_c = k1/rho1/c1;  
    case '2' % Mur tri-couche
        alpha_c = max([k1/rho1/c1 k2/rho2/c2 k3/rho3/c3]);  
    case '3' % Mur avec conductivite en loi de puissance k=b*T^a
        alpha_c = b*max([Tleft, Tmid, Tright])^a/rho1/c1;
    case '4' 
        alpha_c = max([kc/rhoc/cc ki/rhoi/ci]);
end

%% Temps 
Nt     = floor(phi*2*alpha_c*t_final*drm2+10); % nombre de pas de temps
dt     = t_final/Nt;                         % pas de temps
t      = linspace(0,t_final,Nt+1);           % vecteur temps t = [0 ... t_final]

%% Champs de temperature
T      = Tmid*ones(N+1,1);  % declaration du champ de temperature
T(1)   = Tleft;             % temperature au premier point (TEMPERATURE AU CENTRE DU CABLE)
T(N+1) = Tright;            % temperature au dernier point température à l'ext (sur la surface du cable au départ)
Tinit  = T;                 % champ de temperature initial

%% Initialisation des proprietees materiaux
rho   = zeros(N+1,1);
c     = zeros(N+1,1);
theta = zeros(N+1,1);
eta = zeros(N+1,1); %IMPORTANT RADIAL SOUUUUUUUUUUUUUUUUUUUURCE
k     = zeros(N+1,1);
for i = 1:N+1
    switch mur
        case '1' % Mur simple
            rho(i) = rho1;
            c(i)   = c1;
            k(i)   = k1;  
        case '2' % Mur tri-couche
            if r(i)<=r12
                rho(i) = rho1;
                c(i)   = c1;
                k(i)   = k1; 
            elseif r(i)>r12 && r(i)<r23
                rho(i) = rho2;
                c(i)   = c2;
                k(i)   = k2; 
            else%if r(i)>=r23
                rho(i) = rho3;
                c(i)   = c3;
                k(i)   = k3; 
            end  
        case '3' % Mur avec conductivite k en loi de puissance k=b*T^a
            rho(i)   = rho1;
            c(i)     = c1;
            k(i)     = b*T(i)^a;
        case '4' 
            if r(i)<=r12
                rho(i) = rhoc;
                c(i)   = cc;
                k(i)   = kc; 
            else%if r(i)>r12
                rho(i) = rhoi;
                c(i)   = ci;
                k(i)   = ki; 
            end 
    end
   
    theta(i) = 0.5*dt/rho(i)/c(i)*drm2; %theta radial
    %ETA POUR SOURCE
    eta(i) = dt/rho(i)/c(i);
end

%% Declarations pour resolution du systeme
A    = zeros(N+1,1);   % Vecteur contenant les A_i^n
B    = zeros(N+1,1);   % Vecteur contenant les B_i^n
C    = zeros(N+1,1);   % Vecteur contenant les C_i^n
D    = zeros(N+1,1);   % Vecteur contenant les D_i^n
M    = zeros(N-1,N-1); % Matrice contenant M^n

%% Boucle sur le temps
time = 0;
for n = 1:Nt
    % Boucle sur l'espace pour calculer :
    % k+1/2, k-1/2, r+1/2, r-1/2, A_i^n, B_i^n, C_i^n et D_i^n
    % Attention, notre numerotation va de 0 a N.
    % Celle de Matlab va de 1 a N+1.
    % Donc "indice Matlab" = "notre indice" +1
    for i = 2:N
        switch average 
            case 'arithmetic' % Moyennes arithmetiques
                km12  = 0.5*(k(i-1)+k(i));
                kp12  = 0.5*(k(i)+k(i+1));
                rm12  = 0.5*(r(i-1)+r(i));
                rp12  = 0.5*(r(i)+r(i+1));
            case 'harmonic'   % Moyennes harmoniques
                km12  = 2*k(i-1)*k(i)/(k(i-1)+k(i));
                kp12  = 2*k(i)*k(i+1)/(k(i)+k(i+1));
                rm12  = 2*r(i-1)*r(i)/(r(i-1)+r(i));
                rp12  = 2*r(i)*r(i+1)/(r(i)+r(i+1));
            case 'geometric'  % Moyennes geometriques
                km12  = sqrt(k(i-1)*k(i));
                kp12  = sqrt(k(i)*k(i+1));
                rm12  = sqrt(r(i-1)*r(i));
                rp12  = sqrt(r(i)*r(i+1));
        end
        
        % matrice et vecteur en radial
        divqn = (rp12*kp12*(T(i+1)-T(i)  ) ...
                 +rm12*km12*(T(i-1)  -T(i)) );
        A(i) =     - theta(i)* (km12*rm12)/(r(i));        % A_i^n
        B(i) = 1.0 + theta(i)* ((km12*rm12+kp12*rp12))/(r(i)); % B_i^n
        C(i) =     - theta(i)* (kp12*rp12)/(r(i));        % C_i^n
        D(i) = T(i) - theta(i)*(1/(r(i)))*divqn;               % D_i^n
        % Ajout source dans D
        D(i) = D(i) + eta(i)*S;
    end

    % Conditions limites radiale flux null à gauche et flux imposé à droite
    % pour matrice D
    D(2) = 0;
    D(3) = D(3) - A(2)*T(1);
    D(N+1) = -h*T(N+1); %h convection
    % Matrice M^n
    M = diag(B(2:N),0)+diag(C(2:N-1),1)+diag(A(3:N),-1); 
    % Conditions limites radiale flux null à gauche et flux imposé à
    % droite pour matrice M
    M(1,1)= (3/2)*(k(1)/dr);
    M(1,2)= -2*(k(1)/dr);
    M(1,3)= (1/2)*(k(1)/dr);
    M(end, end)= (-3/2)*(k(end)/dr) -h;
    M(end, end-1)= 2*(k(end)/dr);
    M(end, end-2)=(-1/2)*(k(end)/dr);
    % Resolution du systeme lineaire
    Tnum = M\D(2:N); % inconnues ne comprenant pas Tleft et Tright
    T = [T(1); Tnum; T(N+1)]; % concatenation avec Tleft et Tright
    
    % Mise a jour de k si k depend de T, sinon k constant
    if mur == '3'
            k(:) = b.*T(:).^a;
    end

    % Compteurs de fin pour affichages
    time = time + dt;
    if mod(n,n_frequency) == 0
        fprintf('Iteration %d/%d and time %f/%f \n',n,Nt,time,t_final)
    end
    %keyboard
end
%% Plot des resultats

figure(2)
plot(r,Tinit, '-', 'Color','b','LineWidth',4);
hold on
plot(r,T, '-', 'Color','r','LineWidth',4);
set(gca,'FontSize',20)
set(gca,'LineWidth',3)
xlabel('$x$','Interpreter','Latex')
ylabel('$T(x)$','Interpreter','Latex')
grid on
legend('$T(t=0,x)$',['$T(t=',num2str(t_final),',x)$'],'Interpreter','Latex','Fontsize',20,'Location', 'Best')
xline(Rc * R, '--k', 'LineWidth', 2, 'DisplayName', 'Limite cuivre/caoutchouc');



toc