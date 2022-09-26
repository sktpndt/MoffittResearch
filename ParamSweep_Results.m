%% Data

% IMMUNE PARAMETERS %

sigma_delta = [0.1, 0.00, 0.30;
               0.2, 0.30, 0.50;
               0.3, 0.50, 0.60;
               0.4, 0.60, 0.70;
               0.5, 0.80, 0.86;
               0.6, 0.90, 0.94;
               0.7, 1.00, 1.02;
               0.8, 1.06, 1.09;
               ];

% unreliable tuning parameters for pushing the separatrix right/left.
% there's not much of an overlap of parameters that give a good separatrix

sigma_rho = [0.1, 1.0, 1.3;
             0.2, 0.7, 1.1;
             0.3, 0.5, 0.9;
             0.4, 0.4, 0.7;
             0.5, 0.3, 0.5;
             0.6, 0.0, 0.4;
             0.7, 0.0, 0.3;
             0.8, 0.0, 0.2;
             0.9, 0.0, 0.08;
             ];

% increase in rho pushes separatrix to the left 
% increase in sigma pushes separatrix to the left

sigma_mu = [0.1, 0.0021, 0.0028;
            0.2, 0.0027, 0.0061;
            0.3, 0.0033, 0.011;
            0.4, 0.0039, 0.018;
            0.5, 0.0046, 0.025;
            0.6, 0.0053, 0.5; % false limit, but setting mu higher than this takes too long to compute
            0.7, 0.0061, 0.5;
            0.8, 0.0069, 0.5;
            0.9, 0.0078, 0.5;
            1.0, 0.0087, 0.5;
            2.0, 0.019, 0.5;
            3.0, 0.04, 0.5;
            10, 0.13, 0.5;
            20, 0.26, 0.5;
            30, 0.39, 1; % still works all the way up to mu = 1 (1min to compute)
            ];

% increases in sigma push the separatrix up
% increases in mu makes the separatrix flatter

delta_rho = [0.0, 0.00, 0.09;
             0.3, 0.79, 1.24;
             0.4, 0.98, 1.34;
             0.5, 1.14, 1.45;
             0.6, 1.29, 1.55;
             0.7, 1.43, 1.66;
             0.8, 1.57, 1.77;
             0.9, 1.70, 1.87;
             1.0, 1.83, 1.98;
             1.3, 2.20, 2.29;
             1.6, 2.54, 2.61;
             1.9, 2.89, 2.92;
             2.1, 3.11, 3.13;
             2.4, 3.44, 3.45;
             2.7, 3.76, 3.77;
             2.8, 3.87, 3.87;
            ];

% also unreliable as a tuning parameter pair, not much overlap
% DON'T USE DELTA for tuning

delta_mu = [0.0, 0.0032, 0.0220;
            0.1, 0.0030, 0.0110;
            0.2, 0.0027, 0.0070;
            0.3, 0.0024, 0.0043;
            0.4, 0.0022, 0.0029;
            0.5, 0.0019, 0.0020;
            ];

% increase in delta pushes the separatrix to the right
% increases in mu flattens the separatrix

rho_mu = [0.6, 0.00123, 0.00126;
          0.7, 0.0016, 0.0017;
          0.9, 0.0021, 0.0028;
          1.0, 0.0024, 0.0036;
          1.1, 0.0026, 0.0045;
          1.2, 0.0029, 0.0054;
          1.5, 0.0036, 0.0090;
          2.0, 0.0050, 0.0160; 
          3.0, 0.0090, 0.0360;
          10.0, 0.0490, 0.2000; 
          20.0, 0.1280, 0.3000;
          ];

% increase in rho pushes separatrix to the left
% increase in mu flattens the separatrix

% TUMOR PARAMETERS %

alpha_beta = [0.69, 0.001, 0.00103;
              0.9, 0.001, 0.002;
              1, 0.001, 0.0023;
              1.2, 0.001, 0.0025;
              1.3, 0.001, 0.0026;
              1.6, 0.0016, 0.0029;
              1.7, 0.0021, 0.0029;
              2, 0.0029, 0.0031;
              2.3, 0.0033, 0.0033;
              ];

% increasing alpha pushes the separatrix down
% increasing beta pushes separatrix to the left
% beta is a very, very sensitive parameter. Cannot deviate much

alpha_gamma = [0.84, 0.78, 0.94;
               0.9, 0.8, 1.05;
               1, 0.67, 1.23;
               1.2, 0.56, 1.34;
               1.5, 0.83, 1.49;
               2, 1.33, 1.73;
               2.5, 1.9, 1.97;
               ];
% increasing alpha pushes separatrix to the right
% increasing gamma curls the separatrix to the left
% gamma is more forgiving than beta

% alpha is an okay parameter to use, although it has limited range

beta_gamma = [0.0018, 1, 1.64;
              0.002, 0.95, 1.56;
              0.0021, 0.93, 1.52;
              0.0022, 0.9, 1.47;
              0.0023, 0.88, 1.43;
              0.0024, 0.85, 1.38;
              0.0025, 0.82, 1.33;
              0.0027, 0.78, 1.2;
              0.0028, 0.75, 1.13;
              0.0029, 0.78, 1.03;
              ];
% increasing beta shift the separatrix to the left
% increasing gamma shifts the separatrix to the left

% separatrix is very sensitive to changes in beta. Might be a good
% parameter to use for tuning

% separatrix is not very sensitive to gamma, and it has a limited range


% TUMOR V IMMUNE SWEEPS %

beta_rho = [0.001, 1.03, 2.00;
            0.0018, 0.95, 1.45;
            0.0019, 0.94, 1.38;
            0.002, 0.93, 1.32;
            0.0025, 0.88, 1.08;
            0.0029, 0.84, 0.95;
            0.003, 0.83, 0.93;
            0.0033, 0.81, 0.86;
            0.0035, 0.79, 0.82;
            0.004, 0.74, 0.75;
            ];

% increases in beta move separatrix up and to the left (no change in y)
% increases in rho shifts up and to the left (change in y)

beta_mu = [0.0005, 0.0006, 0.002;
           0.001, 0.0012, 0.0023;
           0.002, 0.0022, 0.0032;
           0.0025, 0.0027, 0.0036;
           0.003, 0.0032, 0.004;
           0.0033, 0.0035, 0.0042;
           0.0035, 0.0037, 0.0044;
           0.004, 0.0042, 0.0047;
           0.0045, 0.0046, 0.0051;
           0.005, 0.0051, 0.0054;
           0.0055, 0.0055, 0.0058;
           0.006, 0.0059, 0.0061;
           0.0065, 0.0063, 0.0064;
           0.007, 0.0067, 0.0067;
           0.0074, 0.007, 0.007;
           ];

% beta is consistent
% increase in mu pushes separatrix down and to the right (no change in y)

gamma_rho = [0.6, 1.14, 1.2;
             1, 0.93, 1.32;
             1.5, 0.75, 1;
             1.7, 0.69, 0.81;
             1.8, 0.66, 0.71;
             ];

% increase in gamma shifts up and to the left (change in y)
% rho is consistent

gamma_mu = [1, 0.0022, 0.0032;
            1.5, 0.003, 0.0054;
            2, 0.0042, 0.0081;
            3, 0.0067, 0.015;
            4, 0.0092, 0.023;
            %10, 0.025, 0.5;
            ];

% gamma is conssistent
% mu is consistent

%% immune plots

% Sigma sweeps
figure(1);clf

subplot(1,3,1)
%subplot(1,1,1);
dat = sigma_delta;
tix = linspace(min(dat(:,2)), max(dat(:,3)), 7);
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\sigma", "FontSize",20)
ylabel('\delta', 'FontSize',20)
%%%

subplot(1,3,2)
%subplot(1,1,1);
dat = sigma_rho;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\sigma", "FontSize",20)
ylabel("\rho", "FontSize",20)
%%%
subplot(1,3,3)
%subplot(1,1,1);
dat = sigma_mu;
boxplot(dat(:, 2:3)', dat(:,1));
% plot(dat(:,1)', dat(:,2)')
% hold on;
% plot(dat(:,1)', dat(:,3)')
xlabel("\sigma", "FontSize",20)
ylabel("\mu", "FontSize",20)
%%

t = sgtitle("\sigma Sweeps")
t.FontSize = 25

% Delta sweeps
%%%
figure(2);clf
subplot(1,2,1)
%subplot(1,1,1)
dat = delta_rho;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\delta", "FontSize",20)
ylabel("\rho", "FontSize",20)
%%%
subplot(1,2,2)
%subplot(1,1,1)
dat = delta_mu;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\delta", "FontSize",20)
ylabel("\mu", "FontSize",20)
%%%
t = sgtitle("\delta Sweeps")
t.FontSize = 25

%%%
% Rho sweep
figure(3);clf
subplot(1,1,1)
dat = rho_mu;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\rho", "FontSize",20)
ylabel("\mu", "FontSize",20)
t = sgtitle("\rho Sweep")
t.FontSize = 25
%%%

%% tumor plots

% alpha sweep
%%%
figure(4);clf

subplot(1,2,1)
dat = alpha_beta;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\alpha", "FontSize",20)
ylabel("\beta", "FontSize",20)

%%%

subplot(1,2,2)
%subplot(1,1,1)
dat = alpha_gamma;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\alpha", "FontSize",20)
ylabel("\gamma", "FontSize",20)
t = sgtitle("\alpha Sweep")
t.FontSize = 25

%beta sweep

figure(5);clf
%%%
subplot(1,1,1)
dat = beta_gamma;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\beta", "FontSize",20)
ylabel("\gamma", "FontSize",20)
%%%
t = sgtitle("\beta Sweep")
t.FontSize = 25


%% tumor v immune plots

figure(6);clf

% subplot(1,1,1)
subplot(2,2,1)
dat = beta_rho;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\beta", "FontSize",20)
ylabel("\rho", "FontSize",20)
%%

% subplot(1,1,1)
subplot(2,2,2)
dat = beta_mu;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\beta", "FontSize",20)
ylabel("\mu", "FontSize",20)
%%

% subplot(1,1,1)
subplot(2,2,3)
dat = gamma_rho;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\gamma", "FontSize",20)
ylabel("\rho", "FontSize",20)
%%

% subplot(1,1,1)
subplot(2,2,4)
dat = gamma_mu;
boxplot(dat(:, 2:3)', dat(:,1));
xlabel("\gamma", "FontSize",20)
ylabel("\mu", "FontSize",20)
%%

t = sgtitle("Tumor v Immune Sweeps")
t.FontSize = 25

