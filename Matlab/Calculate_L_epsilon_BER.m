close all;
clear all;

%% calculation parameters
file_version = 2;
ode = 45; % type of ordinary differential equation
filter_type = 'RRC'; %'RRC'; % 'RRC' root raised cosine filter, 'RC' raised cosine filter
Nsym = 6; % RC or RRC filter span in symbol durations
beta = 0.35; % RC or RRC filter roll-off factor
sampsPerSym = 1000; % samples per symbol = Fs/R; must be an integer
R = 1e6; % symbol rate, Hz
Fs = sampsPerSym * R; % sampling rate, usually 1 GHz
Nsymbol = 1000; % number of symbols in the stream
theta_step_size = 1; % angular resolution for the symbol cloud rotation in the constellation diagram, deg
filter_coef = 3; % fc1 = filter_coef*R for demodulation LPF, choose 3 or 10

%% noise between transmitter and receiver and after receiver 
deltaF = 0; % white Gaussian noise in the driving force waveform, N
deltax = 0; % white Gaussian noise in the mechanical response x(t), m
deltarcbackground = 0; % white Gaussian noise at the output of receiver, m

%% mechanical parameters
f0 = 56e6; % mechanical resonant frequency, Hz
fmin = 54e6; % lower bound of the drive frequency range, Hz 
fmax = 58e6; % upper bound of the drive frequency range, Hz
Nfreqoffset = 41; % number of drive frequencies
alpha = -1e14; % Duffing parameter, kg m^(-2) s^(-2)
eta = 1.2e5; % nonlinear damping coefficient, kg m^(-2) s^(-1)
Q0 = 100; % linear quality factor (unitless)
meff = 1.829194458835661e-17; % modal effective mass, kg
Vg = 13; % dc gate voltage, V
r = 3.1e-6/2; % radius of cavity, m     
d = 270e-9;  % cavity depth, m           
epsilon0 = 8.854187817e-12; % permittivity of vacuum, F/m
dC_dz = epsilon0*pi*r^2/d^2; % gradient of gate capacitance, F/m
PdBm = -34; % drive power in dBm

%%
myName = mfilename;

if strcmp(filter_type,'RC') == 1
rctFilt = comm.RaisedCosineTransmitFilter(Shape='Normal',RolloffFactor=beta,FilterSpanInSymbols=Nsym,OutputSamplesPerSymbol=sampsPerSym);
% Normalize to obtain maximum filter tap value of 1
b = coeffs(rctFilt);
rctFilt.Gain = 1/max(b.Numerator);
end

if strcmp(filter_type,'RRC') == 1
    % transmit RRC filter
rctFilt3 = comm.RaisedCosineTransmitFilter(...
  Shape='Square root', ...
  RolloffFactor=beta, ...
  FilterSpanInSymbols=Nsym, ...
  OutputSamplesPerSymbol=sampsPerSym);
  brct = coeffs(rctFilt3);
  rctFilt3.Gain = 1/max(brct.Numerator);

    % receive RRC filter
rcrFilt = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'InputSamplesPerSymbol',  sampsPerSym, ...
  'DecimationFactor',       1);
  brcr = coeffs(rcrFilt);
  rcrFilt.Gain = max(brcr.Numerator);
end

%% Generate symbol components I_{in} and Q_{in}
Iraw = sign(randn(Nsymbol,1));
Qraw = sign(randn(Nsymbol,1));
A = mod(atan2d(Qraw,Iraw),360); % phase angles at baseband
diffA = mod(diff(A),360); % for differential encoding

to = (0: size(Iraw,1)*sampsPerSym - 1) / Fs;
dt = to(2)-to(1);
fn = Fs/2;
fltDelay = Nsym / (2*R); % Filter group delay, since raised cosine filter is linear phase and symmetric.

%% RRC filtering to create \tilde{I}_{in}(t) and \tilde{Q}_{in}(t)
if strcmp(filter_type,'RRC') == 1
yoI = rctFilt3([Iraw; zeros(Nsym/2,1)]);
% Correct for propagation delay by removing filter transients
yoI = yoI(round(fltDelay*Fs+1):end);
dec_yoI = decim_joel(yoI,sampsPerSym,0); % with RC filter, this is the same as Iraw

yoQ = rctFilt3([Qraw; zeros(Nsym/2,1)]); 
% Correct for propagation delay by removing filter transients
yoQ = yoQ(round(fltDelay*Fs+1):end);
dec_yoQ = decim_joel(yoQ,sampsPerSym,0); % with RC filter, this is the same as Qraw
end

%% RC filtering to create \tilde{I}_{in}(t) and \tilde{Q}_{in}(t)
if strcmp(filter_type,'RC') == 1
yoI = rctFilt([Iraw; zeros(Nsym/2,1)]);
% Correct for propagation delay by removing filter transients
yoI = yoI(round(fltDelay*Fs+1):end);
dec_yoI = decim_joel(yoI,sampsPerSym,0); % with RC filter, this is the same as Iraw

yoQ = rctFilt([Qraw; zeros(Nsym/2,1)]); 
% Correct for propagation delay by removing filter transients
yoQ = yoQ(round(fltDelay*Fs+1):end);
dec_yoQ = decim_joel(yoQ,sampsPerSym,0); % with RC filter, this is the same as Qraw
end

%% Lowpass filter for demodulation
fc1 = filter_coef*R; % lowpass filter cutoff frequency 
order= 6;
[z1,p1,k1] = butter(order,fc1/fn,'low');
[sos1,g1] = zp2sos(z1,p1,k1);

EVM_opt_mode1 = nan(1,Nfreqoffset);
EVM_opt_mode2 = nan(1,Nfreqoffset);
EVM_opt_mode3 = nan(1,Nfreqoffset);

Lref_opt = nan(1,Nfreqoffset);

dqpsk_BER = nan(1,Nfreqoffset);
dqpsk_BER_rot = nan(1,Nfreqoffset);

omega0 = 2*pi*f0;
fdrive_pick_array = linspace(fmin,fmax,Nfreqoffset);

if alpha<0
    fdrive_array = 2 * f0 - fdrive_pick_array;
else
    fdrive_array = fdrive_pick_array;
end

alpha = abs(alpha); % alpha>0 is easier on ODE 45
varU = 10^(PdBm/10)*50*1e-3; % V^2_rms
x_init = [0 0];

%% start frequency scan
for w=1:Nfreqoffset
tic
fdrive_pick_array(w)
fdrive = fdrive_array(w);
omega_d = 2*pi*fdrive;

%% Build modulated waveform
x_raw = yoI.*cos(omega_d*to') + yoQ.*sin(omega_d*to');
var_xraw = var(x_raw);
Vpk = sqrt(varU/var_xraw);
Famp = abs(Vg)*2*Vpk*dC_dz; % N, factor of 2 because of full reflection

%% solve equation of motion, set deltaF=0 for noiseless driving force, make sure Famp~=0
[t,x] = eq_motion_JM_QPSK_v2B(x_init,to,omega0,Q0,meff,alpha,eta,Famp,x_raw+deltaF/Famp*randn(size(x_raw,1),1),dt,ode);

%% demodulate \tilde{I}_{out} and \tilde{Q}_{out} from the mechanical response, deltax = 0 for noiseless displacement waveform
I_lp_B = filtfilt(sos1,g1,(x(:,1) + deltax*randn(size(x,1),1)).*cos(omega_d*to'));
Q_lp_B = filtfilt(sos1,g1,(x(:,1) + deltax*randn(size(x,1),1)).*sin(omega_d*to'));

%% RRC filter at the receiver
if strcmp(filter_type,'RRC') == 1
I_lp = rcrFilt([I_lp_B; zeros(Nsym*sampsPerSym/2, 1)]);
% Correct for propagation delay by removing filter transients
I_lp = I_lp(round(fltDelay*Fs+1):end);
I_lp = I_lp + deltarcbackground*randn(size(I_lp,1),1); % deltarcbackground = 0 for noiseless reception

Q_lp = rcrFilt([Q_lp_B; zeros(Nsym*sampsPerSym/2, 1)]);
% Correct for propagation delay by removing filter transients
Q_lp = Q_lp(round(fltDelay*Fs+1):end);
Q_lp = Q_lp + deltarcbackground*randn(size(Q_lp,1),1); % deltarcbackground = 0 for noiseless reception
end

%% Only RC filter
if strcmp(filter_type,'RC') == 1
I_lp = I_lp_B;
I_lp = I_lp + deltarcbackground*randn(size(I_lp,1),1); % deltarcbackground = 0 for noiseless reception

Q_lp = Q_lp_B;
Q_lp = Q_lp + deltarcbackground*randn(size(Q_lp,1),1); % deltarcbackground = 0 for noiseless reception
end

%% Optimize decimation offset of \tilde{I}_{out} and \tilde{Q}_{out}
Lref = nan(1,floor(sampsPerSym/2));
davg_mode1 = nan(1,floor(sampsPerSym/2));
davg_mode2 = nan(1,floor(sampsPerSym/2));
davg_mode3 = nan(1,floor(sampsPerSym/2));
EVM_mode1 = nan(1,floor(sampsPerSym/2));
EVM_mode2 = nan(1,floor(sampsPerSym/2));
EVM_mode3 = nan(1,floor(sampsPerSym/2));

for q=1:1:floor(sampsPerSym/2)
  
decim_offset = q-1;

[dec_I_lp_rot,dec_Q_lp_rot,langle] = rotate_IQ(decim_joel(I_lp,sampsPerSym,decim_offset),decim_joel(Q_lp,sampsPerSym,decim_offset),theta_step_size);

[Lref(q),davg_mode1(q),davg_mode2(q),davg_mode3(q)] = EVM_joel_v6(dec_I_lp_rot,dec_Q_lp_rot);
EVM_mode1(q) = 20*log10(davg_mode1(q)/Lref(q));
EVM_mode2(q) = 20*log10(davg_mode2(q)/Lref(q));
EVM_mode3(q) = 20*log10(davg_mode3(q)/Lref(q));
end

q_opt = find(EVM_mode2==min(EVM_mode2)) -1 ;
EVM_opt_mode2(w) = EVM_mode2(q_opt+1);
EVM_opt_mode1(w) = EVM_mode1(q_opt+1);
EVM_opt_mode3(w) = EVM_mode3(q_opt+1);
Lref_opt(w) = Lref(q_opt+1);

dec_I_lp_opt = decim_joel(I_lp,sampsPerSym,q_opt);
dec_Q_lp_opt = decim_joel(Q_lp,sampsPerSym,q_opt);
[dec_I_lp_rot_opt,dec_Q_lp_rot_opt,langleopt] = rotate_IQ(dec_I_lp_opt,dec_Q_lp_opt,theta_step_size);

%% get symbol error rate
dec_I_lp_opt_norm = nan(1,size(dec_I_lp_opt,2));
dec_Q_lp_opt_norm = nan(1,size(dec_Q_lp_opt,2));

dec_I_lp_rot_opt_norm = nan(1,size(dec_I_lp_rot_opt,2));
dec_Q_lp_rot_opt_norm = nan(1,size(dec_Q_lp_rot_opt,2));

for n=1:size(dec_I_lp_opt,2)
      dec_I_lp_opt_norm(n) = sign(dec_I_lp_opt(1,n));
      dec_Q_lp_opt_norm(n) = sign(dec_Q_lp_opt(1,n));
      dec_I_lp_rot_opt_norm(n) = sign(dec_I_lp_rot_opt(1,n));
      dec_Q_lp_rot_opt_norm(n) = sign(dec_Q_lp_rot_opt(1,n));
end

% differential QPSK unrotated
C = mod(atan2d(dec_Q_lp_opt_norm,dec_I_lp_opt_norm),360);
diffC = mod(diff(C'),360);

ladiff = mod(diffA-diffC,180);
ladiff = ladiff(2:end,1);
zerodiff = ladiff(ladiff == 0);
dqpsk_BER(w) = (1 - size(zerodiff,1)/size(ladiff,1))*100; % percent

% differential QPSK rotated
Crot = mod(atan2d(dec_Q_lp_rot_opt_norm,dec_I_lp_rot_opt_norm),360);
diffCrot = mod(diff(Crot'),360);

ladiffrot = mod(diffA-diffCrot,180);
ladiffrot = ladiffrot(2:end,1);
zerodiffrot = ladiffrot(ladiffrot == 0);
dqpsk_BER_rot(w) = (1 - size(zerodiffrot,1)/size(ladiffrot,1))*100; % percent

toc
end

%% Below, \tilde{I}_{out}(t), \tilde{Q}_{out}(t), {I}_{out}, Q_{out}, and L are scaled by a factor of 2 because they are obtained from downconverting the response to baseband

hfig1 = figure;

subplot(2,1,1)
plot(2*dec_I_lp_opt,2*dec_Q_lp_opt,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3)
xlabel('I_{out} (m)');
ylabel('Q_{out} (m)');

subplot(2,1,2)
plot(2*dec_I_lp_rot_opt,2*dec_Q_lp_rot_opt,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3)
xlabel('Rotated I_{out} (m)');
ylabel('Rotated Q_{out} (m)');

set(gcf, 'Position',  [1200, 0, 350, 1200]);

hfig2 = figure;

subplot(4,1,1)
line([0 to(end)],[0 0],'LineStyle','--','color','k')
hold on
plot(to,yoI,'linewidth',1,'color',[1 0 0]);
hold on
stem(linspace(0,size(Iraw,1)-1,size(Iraw,1))/R,Iraw, 'kx');
xlabel('t (s)');
ylabel('I_{in}');
xlim([0 to(end)])
ylim([-1.1*max(max(abs(yoI),max(abs(yoQ)))),1.1*max(max(abs(yoI),max(abs(yoQ))))])
box on;

subplot(4,1,2)
line([0 to(end)],[0 0],'LineStyle','--','color','k')
hold on
plot(to,yoQ,'linewidth',1,'color',[0 0 1]);
hold on
stem(linspace(0,size(Qraw,1)-1,size(Qraw,1))/R,Qraw, 'kx');
xlabel('t (s)');
ylabel('Q_{in}');
xlim([0 to(end)])
ylim([-1.1*max(max(abs(yoI),max(abs(yoQ)))),1.1*max(max(abs(yoI),max(abs(yoQ))))])
box on;

subplot(4,1,3)
line([0 to(end)],[0 0],'LineStyle','--','color','k')
hold on
plot(to,2*I_lp,'linewidth',2,'color',[1 0.6 0.6]);
hold on
stem(decim_joel(to',sampsPerSym,q_opt),2*dec_I_lp_opt, 'kx');
xlabel('t (s)');
ylabel('I_{out} (m)');
xlim([0 to(end)])
ylim([-1.1*max(max(2*abs(I_lp)),max(2*abs(Q_lp))) 1.1*max(max(2*abs(I_lp)),max(abs(2*Q_lp)))])
box on;

subplot(4,1,4)
line([0 to(end)],[0 0],'LineStyle','--','color','k')
hold on
plot(to,2*Q_lp,'linewidth',2,'color',[0.6 0.6 1]);
hold on
stem(decim_joel(to',sampsPerSym,q_opt),2*dec_Q_lp_opt, 'kx');
xlabel('t (s)');
ylabel('Q_{out} (m)');
xlim([0 to(end)])
ylim([-1.1*max(max(2*abs(I_lp)),max(2*abs(Q_lp))) 1.1*max(max(2*abs(I_lp)),max(2*abs(Q_lp)))])
box on;

set(gcf, 'Position',  [0, 0, 1200, 800]);

hfig3 = figure;
subplot(3,1,1)
plot(fdrive_pick_array,2*Lref_opt,'-o','linewidth',1,'color','b');
xlabel('f_d (Hz)')
ylabel('L (m)')
xlim([fmin fmax]);

subplot(3,1,2)
plot(fdrive_pick_array,-EVM_opt_mode2,'-o','linewidth',1,'color','r');
xlabel('f_d (Hz)')
ylabel('1/EVM (dB)')
xlim([fmin fmax]);

subplot(3,1,3)
plot(fdrive_pick_array,dqpsk_BER/100,'-o','linewidth',1,'color','b');
hold on
plot(fdrive_pick_array,dqpsk_BER_rot/100,'-o','linewidth',1,'color','r');
xlabel('f_d (Hz)')
ylabel('BER')
xlim([fmin fmax]);

set(gcf, 'Position',  [800, 0, 400, 800]);
evm_data_path = 'none';
length_data_path = 'none';

result = transpose(cat(1,fdrive_pick_array,2*Lref_opt,EVM_opt_mode1,EVM_opt_mode2,EVM_opt_mode3,dqpsk_BER,dqpsk_BER_rot));
metadata = {strcat('PdBm=',num2str(PdBm),'dBm');strcat('Rsym=',num2str(R/1e6),'MHz');strcat('filter type=',filter_type);strcat('beta=',num2str(beta));strcat('Nsym=',num2str(Nsym));strcat('sampsPerSym=',num2str(sampsPerSym));strcat('f0=',num2str(f0),'Hz');strcat('NumberOfSymbols=',num2str(Nsymbol));strcat('Fs=',num2str(Fs));strcat('filter_coef=',num2str(filter_coef));strcat('deltaF=',num2str(deltaF));strcat('deltax=',num2str(deltax));strcat('deltarcbackground=',num2str(deltarcbackground));strcat('theta_step_size=',num2str(theta_step_size),'deg');strcat('fmin=',num2str(fmin),'Hz');strcat('fmax=',num2str(fmax),'Hz');strcat('Nfreqoffset=',num2str(Nfreqoffset));strcat('alpha=',num2str(alpha));strcat('eta=',num2str(eta));strcat('Q0=',num2str(Q0));strcat('ODE',num2str(ode));strcat('sim version: ',myName)};

lenom = strcat('.\P',num2str(PdBm),'dBm_Rp',num2str(R/1e6),'MHz_',strcat('v',num2str(file_version)));
save(strcat(lenom,'.dat'),'result','-ascii');
writelines(string(metadata),strcat(lenom,'_meta.txt'));
print(gcf,strcat('.\R',num2str(R/1e6),'MHz_Pm',num2str(PdBm),'dBm_',num2str(file_version),'.png'),'-dpng','-r1200')
