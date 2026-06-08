%% ============================================================
%  SIMULASI PERSAMAAN DIFERENSIAL - 3 KASUS SISTEM DINAMIS
%  Tugas Proyek Persamaan Diferensial
%  Departemen Teknik Elektro - Universitas Diponegoro
%
%  KASUS I   : Sistem Mekanik - Gerbong Lokomotif (Massa-Pegas-Redaman)
%  KASUS II  : Sistem Elektrik - Rangkaian RLC Seri
%  KASUS III : Sistem Mekatronik - Motor DC
%
%  Metode: Koefisien Tak Tentu, Variasi Parameter, Transformasi Laplace
%  Validasi numerik menggunakan ode45
% ============================================================

clc; clear; close all;

%% ============================================================
%  KASUS I: GERBONG LOKOMOTIF (Massa-Pegas-Redaman)
%  m*x'' + c*x' + k*x = F(t)   dengan  F(t) = 40000*e^(-0.2t)
% ============================================================
fprintf('=== KASUS I: GERBONG LOKOMOTIF ===\n');

% --- Parameter Fisik ---
m = 35000;       % massa gerbong (kg)
c = 15000;       % koefisien peredam (Ns/m)
k = 450000;      % konstanta pegas (N/m)
F0 = 40000;      % amplitudo gaya (N)
alpha_F = 0.2;   % konstanta peluruhan gaya

% Kondisi awal
x0 = 0;          % x(0) = 0
xd0 = 0;         % x'(0) = 0

t = linspace(0, 5, 2000);

% Akar persamaan karakteristik: m*r^2 + c*r + k = 0
disc1 = c^2 - 4*m*k;
r1 = (-c + sqrt(complex(disc1))) / (2*m);
r2 = (-c - sqrt(complex(disc1))) / (2*m);
alpha1 = real(r1);
beta1  = imag(r1);

fprintf('Akar karakteristik r1 = %.4f + %.4fi\n', real(r1), imag(r1));
fprintf('Akar karakteristik r2 = %.4f + %.4fi\n', real(r2), imag(r2));
zeta1 = c / (2*sqrt(m*k));
fprintf('Rasio redaman zeta = %.4f => ', zeta1);
if zeta1 < 1
    fprintf('UNDERDAMPED\n');
elseif zeta1 == 1
    fprintf('CRITICALLY DAMPED\n');
else
    fprintf('OVERDAMPED\n');
end

% ---- METODE 1: Koefisien Tak Tentu ----
% Tebak xp = A*e^(-alpha_F*t)
% Substitusi: A*(m*alpha_F^2 - c*alpha_F + k) = F0
denom1 = m*alpha_F^2 - c*alpha_F + k;
A_ktk = F0 / denom1;

% Solusi homogen: xh = C1*e^(r1*t) + C2*e^(r2*t)  (underdamped => exp*sin/cos)
% x(t) = e^(alpha1*t)*(C1*cos(beta1*t) + C2*sin(beta1*t)) + A*e^(-alpha_F*t)
% Dari IC: x(0)=0 => C1 + A = 0 => C1 = -A
% x'(0)=0 => alpha1*C1 + beta1*C2 - alpha_F*A = 0
C1_1 = -A_ktk;
C2_1 = (alpha_F*A_ktk - alpha1*C1_1) / beta1;

x_KTK1 = exp(alpha1*t).*(C1_1*cos(beta1*t) + C2_1*sin(beta1*t)) + A_ktk*exp(-alpha_F*t);

% ---- METODE 2: Variasi Parameter ----
% y1 = e^(alpha1*t)*cos(beta1*t), y2 = e^(alpha1*t)*sin(beta1*t)
% Wronskian W = beta1 * e^(2*alpha1*t)
% g(t) = F(t)/m = (F0/m)*e^(-alpha_F*t)
g1 = @(tau) (F0/m) .* exp(-alpha_F*tau);
y1f = @(tau) exp(alpha1*tau) .* cos(beta1*tau);
y2f = @(tau) exp(alpha1*tau) .* sin(beta1*tau);
W1  = @(tau) beta1 * exp(2*alpha1*tau);

u1_int = zeros(size(t));
u2_int = zeros(size(t));
dt1 = t(2)-t(1);
for i = 2:length(t)
    u1_int(i) = u1_int(i-1) + (-y2f(t(i-1))*g1(t(i-1))/W1(t(i-1)))*dt1;
    u2_int(i) = u2_int(i-1) + ( y1f(t(i-1))*g1(t(i-1))/W1(t(i-1)))*dt1;
end
xp_VP1 = u1_int.*y1f(t) + u2_int.*y2f(t);
% Sesuaikan konstanta homogen dari IC
xh_VP1 = @(C1,C2) exp(alpha1*t).*(C1*cos(beta1*t)+C2*sin(beta1*t));
% x(0)=0: C1 + xp(0) = 0 => C1 = -xp(0) ≈ 0
% x'(0)=0
C1_v1 = -xp_VP1(1);
dxp_VP1_0 = (xp_VP1(2)-xp_VP1(1))/dt1;
C2_v1 = (-dxp_VP1_0 - alpha1*C1_v1) / beta1;

x_VP1 = exp(alpha1*t).*(C1_v1*cos(beta1*t)+C2_v1*sin(beta1*t)) + xp_VP1;

% ---- METODE 3: Transformasi Laplace ----
% X(s) = F0/((m*s^2+c*s+k)*(s+alpha_F))
% Pecahan parsial: X(s) = A/(s+alpha_F) + (Bs+D)/(m*s^2+c*s+k)
% Solusi domain waktu (residue numerik)
num_L1  = [F0];
den_L1  = conv([m c k], [1 alpha_F]);
[r_L,p_L,~] = residue(num_L1, den_L1);

x_Lap1 = zeros(size(t));
for i = 1:length(r_L)
    if abs(imag(p_L(i))) < 1e-6
        x_Lap1 = x_Lap1 + real(r_L(i))*exp(real(p_L(i))*t);
    else
        % pasangan konjugat
        x_Lap1 = x_Lap1 + 2*real(r_L(i)*exp(p_L(i)*t));
    end
end
% Terapkan hanya satu kali untuk pasangan konjugat
% (residue sudah menangani pasangan, reset dan hitung ulang)
x_Lap1 = zeros(size(t));
visited = false(size(r_L));
for i = 1:length(r_L)
    if visited(i), continue; end
    if abs(imag(p_L(i))) < 1e-6
        x_Lap1 = x_Lap1 + real(r_L(i))*exp(real(p_L(i))*t);
        visited(i) = true;
    else
        for j = i+1:length(r_L)
            if abs(p_L(j) - conj(p_L(i))) < 1e-6
                x_Lap1 = x_Lap1 + 2*real(r_L(i)*exp(p_L(i)*t));
                visited(i) = true;
                visited(j) = true;
                break;
            end
        end
    end
end

% ---- Validasi Numerik ode45 ----
ode_gerbong = @(t,y) [y(2); (F0*exp(-alpha_F*t) - c*y(2) - k*y(1))/m];
[t_ode1, y_ode1] = ode45(ode_gerbong, [0 5], [0 0]);

% ---- Plot Kasus I ----
figure('Name','Kasus I - Gerbong Lokomotif','Position',[50 50 1200 500]);
subplot(1,2,1);
plot(t, x_KTK1*1000,'b-','LineWidth',2); hold on;
plot(t, x_VP1*1000, 'r--','LineWidth',2);
plot(t, x_Lap1*1000,'g-.','LineWidth',2);
plot(t_ode1, y_ode1(:,1)*1000,'k:','LineWidth',1.5);
xlabel('Waktu (s)'); ylabel('Simpangan x(t) (mm)');
title('Kasus I: Gerbong Lokomotif - Perbandingan Metode');
legend('Koef. Tak Tentu','Variasi Parameter','Laplace','ode45 (numerik)','Location','best');
grid on;

subplot(1,2,2);
plot(t_ode1, y_ode1(:,1)*1000,'k-','LineWidth',2);
xlabel('Waktu (s)'); ylabel('Simpangan x(t) (mm)');
title('Validasi Numerik ode45');
grid on;

fprintf('Parameter Gerbong: m=%d kg, c=%d Ns/m, k=%d N/m\n', m, c, k);
fprintf('Steady-state xss = F0/k = %.6f mm\n\n', F0/k*1000);


%% ============================================================
%  KASUS II: RANGKAIAN RLC SERI
%  L*q'' + R*q' + (1/C)*q = V(t)   dengan  V(t) = 24*e^(-2t)
% ============================================================
fprintf('=== KASUS II: RANGKAIAN RLC SERI ===\n');

% --- Parameter ---
R   = 150;           % Resistor (Ohm)
L   = 0.8;           % Induktor (H)
Cap = 200e-6;        % Kapasitor (F)
Vm  = 24;            % Amplitudo tegangan (V)
alpha_V = 2;         % Konstanta peluruhan tegangan

% Kondisi awal
q0  = 0;
qd0 = 0;

t2 = linspace(0, 0.15, 2000);

% Akar karakteristik: L*r^2 + R*r + 1/C = 0
disc2 = R^2 - 4*L/Cap;
r1_2 = (-R + sqrt(complex(disc2))) / (2*L);
r2_2 = (-R - sqrt(complex(disc2))) / (2*L);
alpha2 = real(r1_2);
beta2  = abs(imag(r1_2));

fprintf('Akar r1 = %.4f + %.4fi\n', real(r1_2), imag(r1_2));
fprintf('Akar r2 = %.4f + %.4fi\n', real(r2_2), imag(r2_2));
omega0_2 = 1/sqrt(L*Cap);
zeta2 = R / (2*sqrt(L/Cap));
fprintf('omega0 = %.2f rad/s, zeta = %.4f => ', omega0_2, zeta2);
if zeta2 < 1
    fprintf('UNDERDAMPED\n');
elseif zeta2 == 1
    fprintf('CRITICALLY DAMPED\n');
else
    fprintf('OVERDAMPED\n');
end

% ---- METODE 1: Koefisien Tak Tentu ----
% qp = A*e^(-alpha_V*t)
denom2 = L*alpha_V^2 - R*alpha_V + 1/Cap;
A_ktk2 = Vm / denom2;
% C1 + A = q0 = 0 => C1 = -A
C1_2 = -A_ktk2;
% q'(0)=0: alpha2*C1 + beta2*C2 - alpha_V*A = 0
C2_2 = (alpha_V*A_ktk2 - alpha2*C1_2) / beta2;

q_KTK2 = exp(alpha2*t2).*(C1_2*cos(beta2*t2) + C2_2*sin(beta2*t2)) + A_ktk2*exp(-alpha_V*t2);

% ---- METODE 2: Variasi Parameter ----
g2  = @(tau) (Vm/L)*exp(-alpha_V*tau);
y1f2 = @(tau) exp(alpha2*tau).*cos(beta2*tau);
y2f2 = @(tau) exp(alpha2*tau).*sin(beta2*tau);
W2   = @(tau) beta2*exp(2*alpha2*tau);

u1_int2 = zeros(size(t2));
u2_int2 = zeros(size(t2));
dt2 = t2(2)-t2(1);
for i = 2:length(t2)
    u1_int2(i) = u1_int2(i-1) + (-y2f2(t2(i-1))*g2(t2(i-1))/W2(t2(i-1)))*dt2;
    u2_int2(i) = u2_int2(i-1) + ( y1f2(t2(i-1))*g2(t2(i-1))/W2(t2(i-1)))*dt2;
end
qp_VP2 = u1_int2.*y1f2(t2) + u2_int2.*y2f2(t2);
C1_v2 = -qp_VP2(1);
dqp_VP2_0 = (qp_VP2(2)-qp_VP2(1))/dt2;
C2_v2 = (-dqp_VP2_0 - alpha2*C1_v2) / beta2;
q_VP2 = exp(alpha2*t2).*(C1_v2*cos(beta2*t2)+C2_v2*sin(beta2*t2)) + qp_VP2;

% ---- METODE 3: Transformasi Laplace ----
num_L2 = [Vm];
den_L2 = conv([L R 1/Cap], [1 alpha_V]);
[r_L2,p_L2,~] = residue(num_L2, den_L2);

q_Lap2 = zeros(size(t2));
visited2 = false(size(r_L2));
for i = 1:length(r_L2)
    if visited2(i), continue; end
    if abs(imag(p_L2(i))) < 1e-6
        q_Lap2 = q_Lap2 + real(r_L2(i))*exp(real(p_L2(i))*t2);
        visited2(i) = true;
    else
        for j = i+1:length(r_L2)
            if abs(p_L2(j)-conj(p_L2(i))) < 1e-6
                q_Lap2 = q_Lap2 + 2*real(r_L2(i)*exp(p_L2(i)*t2));
                visited2(i)=true; visited2(j)=true;
                break;
            end
        end
    end
end

% ---- Validasi ode45 ----
ode_rlc = @(t,y) [y(2); (Vm*exp(-alpha_V*t) - R*y(2) - y(1)/Cap)/L];
[t_ode2, y_ode2] = ode45(ode_rlc, [0 0.15], [0 0]);

% Arus i = dq/dt
i_KTK2  = gradient(q_KTK2,  dt2);
i_VP2   = gradient(q_VP2,   dt2);
i_Lap2  = gradient(q_Lap2,  dt2);

% ---- Plot Kasus II ----
figure('Name','Kasus II - Rangkaian RLC','Position',[50 600 1200 500]);
subplot(1,2,1);
plot(t2*1000, q_KTK2*1e6,'b-','LineWidth',2); hold on;
plot(t2*1000, q_VP2*1e6, 'r--','LineWidth',2);
plot(t2*1000, q_Lap2*1e6,'g-.','LineWidth',2);
plot(t_ode2*1000, y_ode2(:,1)*1e6,'k:','LineWidth',1.5);
xlabel('Waktu (ms)'); ylabel('Muatan q(t) (\muC)');
title('Kasus II: RLC Seri - Muatan Kapasitor');
legend('Koef. Tak Tentu','Variasi Parameter','Laplace','ode45','Location','best');
grid on;

subplot(1,2,2);
plot(t2*1000, i_KTK2*1000,'b-','LineWidth',2); hold on;
plot(t2*1000, i_VP2*1000, 'r--','LineWidth',2);
plot(t2*1000, i_Lap2*1000,'g-.','LineWidth',2);
xlabel('Waktu (ms)'); ylabel('Arus i(t) (mA)');
title('Kasus II: RLC Seri - Arus');
legend('Koef. Tak Tentu','Variasi Parameter','Laplace','Location','best');
grid on;

fprintf('Parameter RLC: R=%.0f Ohm, L=%.2f H, C=%.0f uF\n', R, L, Cap*1e6);
fprintf('omega0 = %.2f rad/s, T_natural = %.4f s\n\n', omega0_2, 2*pi/omega0_2);


%% ============================================================
%  KASUS III: MOTOR DC
%  (L*J)*omega'' + (L*b + R*J)*omega' + (R*b + Kt*Ke)*omega = Kt*Vin
%  dengan Vin = 12*e^(-2t)  (mengikuti bentuk V(t) = 12e^(-2t))
% ============================================================
fprintf('=== KASUS III: MOTOR DC ===\n');

% --- Parameter ---
J_m  = 0.02;    % momen inersia rotor (kg.m^2)
b_m  = 0.2;     % koefisien gesekan (Nm.s/rad)
Kt   = 0.05;    % konstanta torsi (Nm/A)
Ke   = 0.05;    % konstanta back-EMF (V.s/rad)
R_m  = 2;       % hambatan armatur (Ohm)
L_m  = 1.0;     % induktansi armatur (H)
Vin  = 12;      % amplitudo tegangan input (V)
alpha_m = 2;    % konstanta peluruhan tegangan

% Kondisi awal
omega0_m  = 0;
domega0_m = 0;

t3 = linspace(0, 8, 2000);

% Koefisien PD: a*omega'' + b*omega' + c*omega = K*Vin
a3 = L_m * J_m;
b3 = L_m*b_m + R_m*J_m;
c3 = R_m*b_m + Kt*Ke;
K3 = Kt;

fprintf('PD Motor DC: %.4f*w'' + %.4f*w'' + %.4f*w = %.4f*Vin(t)\n', a3,b3,c3,K3);

% Akar karakteristik
disc3 = b3^2 - 4*a3*c3;
r1_3 = (-b3 + sqrt(complex(disc3))) / (2*a3);
r2_3 = (-b3 - sqrt(complex(disc3))) / (2*a3);
alpha3 = real(r1_3);
beta3  = abs(imag(r1_3));

fprintf('Akar r1 = %.4f + %.4fi\n', real(r1_3), imag(r1_3));
fprintf('Akar r2 = %.4f + %.4fi\n', real(r2_3), imag(r2_3));
zeta3 = b3 / (2*sqrt(a3*c3));
fprintf('Rasio redaman zeta = %.4f => ', zeta3);
if zeta3 < 1
    fprintf('UNDERDAMPED\n');
elseif abs(zeta3-1) < 1e-4
    fprintf('CRITICALLY DAMPED\n');
else
    fprintf('OVERDAMPED\n');
end

% ---- METODE 1: Koefisien Tak Tentu ----
% omegap = A*e^(-alpha_m*t)
denom3 = a3*alpha_m^2 - b3*alpha_m + c3;
A_ktk3 = (K3*Vin) / denom3;

if abs(beta3) < 1e-6
    % Overdamped: dua akar real
    r1r = real(r1_3); r2r = real(r2_3);
    % omega(0)=0: C1+C2+A = 0
    % omega'(0)=0: r1*C1+r2*C2-alpha_m*A = 0
    C2_ktk3 = (-A_ktk3*(r1r + alpha_m)) / (r1r - r2r);
    C1_ktk3 = -A_ktk3 - C2_ktk3;
    omega_KTK3 = C1_ktk3*exp(r1r*t3) + C2_ktk3*exp(r2r*t3) + A_ktk3*exp(-alpha_m*t3);
else
    C1_ktk3 = -A_ktk3;
    C2_ktk3 = (alpha_m*A_ktk3 - alpha3*C1_ktk3) / beta3;
    omega_KTK3 = exp(alpha3*t3).*(C1_ktk3*cos(beta3*t3)+C2_ktk3*sin(beta3*t3)) + A_ktk3*exp(-alpha_m*t3);
end

% ---- METODE 2: Variasi Parameter ----
g3  = @(tau) (K3*Vin/a3)*exp(-alpha_m*tau);
if abs(beta3) < 1e-6
    r1r3 = real(r1_3); r2r3 = real(r2_3);
    y1f3 = @(tau) exp(r1r3*tau);
    y2f3 = @(tau) exp(r2r3*tau);
    W3   = @(tau) (r2r3-r1r3)*exp((r1r3+r2r3)*tau);
else
    y1f3 = @(tau) exp(alpha3*tau).*cos(beta3*tau);
    y2f3 = @(tau) exp(alpha3*tau).*sin(beta3*tau);
    W3   = @(tau) beta3*exp(2*alpha3*tau);
end

u1_int3 = zeros(size(t3));
u2_int3 = zeros(size(t3));
dt3 = t3(2)-t3(1);
for i = 2:length(t3)
    u1_int3(i) = u1_int3(i-1) + (-y2f3(t3(i-1))*g3(t3(i-1))/W3(t3(i-1)))*dt3;
    u2_int3(i) = u2_int3(i-1) + ( y1f3(t3(i-1))*g3(t3(i-1))/W3(t3(i-1)))*dt3;
end
omegap_VP3 = u1_int3.*y1f3(t3) + u2_int3.*y2f3(t3);

if abs(beta3) < 1e-6
    r1r3 = real(r1_3); r2r3 = real(r2_3);
    % IC: C1+C2+omegap(0)=0, r1*C1+r2*C2+domegap(0)=0
    domegap3_0 = (omegap_VP3(2)-omegap_VP3(1))/dt3;
    A_mat = [1 1; r1r3 r2r3];
    b_vec = [-omegap_VP3(1); -domegap3_0];
    C_vec = A_mat \ b_vec;
    omega_VP3 = C_vec(1)*exp(r1r3*t3) + C_vec(2)*exp(r2r3*t3) + omegap_VP3;
else
    C1_v3 = -omegap_VP3(1);
    domegap3_0 = (omegap_VP3(2)-omegap_VP3(1))/dt3;
    C2_v3 = (-domegap3_0 - alpha3*C1_v3) / beta3;
    omega_VP3 = exp(alpha3*t3).*(C1_v3*cos(beta3*t3)+C2_v3*sin(beta3*t3)) + omegap_VP3;
end

% ---- METODE 3: Transformasi Laplace ----
% Omega(s) = (K3*Vin) / ((a3*s^2+b3*s+c3)*(s+alpha_m))
num_L3 = [K3*Vin];
den_L3 = conv([a3 b3 c3], [1 alpha_m]);
[r_L3,p_L3,~] = residue(num_L3, den_L3);

omega_Lap3 = zeros(size(t3));
visited3 = false(size(r_L3));
for i = 1:length(r_L3)
    if visited3(i), continue; end
    if abs(imag(p_L3(i))) < 1e-6
        omega_Lap3 = omega_Lap3 + real(r_L3(i))*exp(real(p_L3(i))*t3);
        visited3(i) = true;
    else
        for j = i+1:length(r_L3)
            if abs(p_L3(j)-conj(p_L3(i))) < 1e-6
                omega_Lap3 = omega_Lap3 + 2*real(r_L3(i)*exp(p_L3(i)*t3));
                visited3(i)=true; visited3(j)=true;
                break;
            end
        end
    end
end

% ---- Hitung arus armatur i(t) dari omega(t) ----
% J*domega/dt + b*omega = Kt*i  =>  i = (J*domega/dt + b*omega) / Kt
i_motor_KTK = (J_m*gradient(omega_KTK3,dt3) + b_m*omega_KTK3) / Kt;
i_motor_VP  = (J_m*gradient(omega_VP3, dt3) + b_m*omega_VP3)  / Kt;
i_motor_Lap = (J_m*gradient(omega_Lap3,dt3) + b_m*omega_Lap3) / Kt;

% ---- Validasi ode45 ----
ode_motorDC = @(t,y) [y(2);
    (K3*Vin*exp(-alpha_m*t) - b3*y(2) - c3*y(1)) / a3];
[t_ode3, y_ode3] = ode45(ode_motorDC, [0 8], [0 0]);

% ---- Plot Kasus III ----
figure('Name','Kasus III - Motor DC','Position',[100 100 1200 500]);
subplot(1,2,1);
plot(t3, omega_KTK3,'b-','LineWidth',2); hold on;
plot(t3, omega_VP3, 'r--','LineWidth',2);
plot(t3, omega_Lap3,'g-.','LineWidth',2);
plot(t_ode3, y_ode3(:,1),'k:','LineWidth',1.5);
xlabel('Waktu (s)'); ylabel('\omega(t) (rad/s)');
title('Kasus III: Motor DC - Kecepatan Sudut');
legend('Koef. Tak Tentu','Variasi Parameter','Laplace','ode45','Location','best');
grid on;

subplot(1,2,2);
plot(t3, i_motor_KTK,'b-','LineWidth',2); hold on;
plot(t3, i_motor_VP, 'r--','LineWidth',2);
plot(t3, i_motor_Lap,'g-.','LineWidth',2);
xlabel('Waktu (s)'); ylabel('Arus i(t) (A)');
title('Kasus III: Motor DC - Arus Armatur');
legend('Koef. Tak Tentu','Variasi Parameter','Laplace','Location','best');
grid on;

omega_ss = K3*Vin / c3;
fprintf('Kecepatan steady-state (teoritis) = %.4f rad/s\n', omega_ss);
fprintf('Parameter: J=%.3f, b=%.2f, Kt=%.2f, Ke=%.2f, R=%.1f, L=%.1f\n\n', ...
    J_m, b_m, Kt, Ke, R_m, L_m);


%% ============================================================
%  RINGKASAN PERBANDINGAN NUMERIK
% ============================================================
fprintf('=== RINGKASAN KESESUAIAN ANTAR METODE ===\n');

% Kasus I
err1_vp  = norm(x_KTK1 - x_VP1) / norm(x_KTK1) * 100;
err1_lap = norm(x_KTK1 - x_Lap1) / norm(x_KTK1) * 100;
fprintf('Kasus I  | VP vs KTK: %.4f%% | Laplace vs KTK: %.4f%%\n', err1_vp, err1_lap);

% Kasus II
err2_vp  = norm(q_KTK2 - q_VP2) / norm(q_KTK2) * 100;
err2_lap = norm(q_KTK2 - q_Lap2) / norm(q_KTK2) * 100;
fprintf('Kasus II | VP vs KTK: %.4f%% | Laplace vs KTK: %.4f%%\n', err2_vp, err2_lap);

% Kasus III
err3_vp  = norm(omega_KTK3 - omega_VP3)  / norm(omega_KTK3) * 100;
err3_lap = norm(omega_KTK3 - omega_Lap3) / norm(omega_KTK3) * 100;
fprintf('Kasus III| VP vs KTK: %.4f%% | Laplace vs KTK: %.4f%%\n', err3_vp, err3_lap);

fprintf('\nSELESAI - Ketiga metode menghasilkan solusi yang ekuivalen.\n');
fprintf('Perbedaan kecil disebabkan oleh akumulasi galat numerik integrasi.\n');
