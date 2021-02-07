% Kyle Acheson - 07/09/2020
% IAM rotationally averaged 1D scattering pattern calculator
% Can be used for X-ray scattering or electron diffraction
% Works for static of time-dependent signals - will detect automatically
% Uses Adams f_functions code for form factors - H,C,N,F,S,I,Xe available.
clear all 
%close all 
tstart = tic;

% Read in geometries in Angstrom
%load('wigner_dist_n100_293k.mat'); % (Natom,Natom,Ntrj,Nts)
%Q = load('tol.dat');
load('extended_trajs_213.mat')
Q = geometries; 
[x,y,Ntraj,Nts] = size(Q);

% Setup Options
FLAGelec = 1; % 0 for X-ray, 1 for UED
tmax = 1000; % max time in fs
exfrac = 10; % excitation fraction in percentage units
%atmnum = [6 6 6 6 6 6 1 1 1 1 1 6 1 1 1 ]; % atoms
atmnum = [6 16 16];
au2ang = 0.52917721092d0;
ang2au = 1/au2ang;
kin = 1.8751e+03*au2ang; % wave vector - cs2 UED
%kin   = 12.d0*au2ang; % incident wave vector au - cs2 Xray
theta_lim = [0 pi]; % range of theta angles
Nq = 100000; % number of points over theta range 
Q = Q*ang2au; % convert geometries to a.u.
tt = linspace(0,tmax,Nts); % set time vector for plotting

theta_min = theta_lim(1);
theta_max = theta_lim(2);
qmin =2.d0*kin*sin(0.5d0*theta_min);
qmax = 2.d0*kin*sin(0.5d0*theta_max);
q = linspace(qmin,qmax,Nq); % linear in q
qAng = q*ang2au; % q in inv Ang
qAng = qAng(qAng<8*pi); % set the limit of q 
Nq = length(qAng); % redefine
q = q(1:Nq);

[FF,fq]  = get_scattering_factors(qAng,atmnum,FLAGelec); % NB: call using q in Angstroms^(-1)
Iat      = sum(fq.^2); % atomic scattering term 

Wiam_tot = zeros(Nq,Nts);
for ts=1:Nts % loop over time steps
    Wiam = zeros(Ntraj,Nq); % init scattering matrix for N trajs
    D = Q(:,:,:,ts); % geometry at each time-step
    for traj=1:Ntraj % loop over trajs
        Imol = zeros(1,Nq); % init for each traj at each time step
        sinQ = zeros(1,Nq); 
        for a=1:x
            for b=a+1:x
                DQ = q(1:Nq)*norm(D(a,1:3,traj)-D(b,1:3,traj)); %qRij
                sinQ(1:Nq) = sin(DQ(1:Nq))./DQ(1:Nq); % Sin(qRij)/qRij
                ind = find(abs(DQ)<1.d-9); % check q=0 limit
                sinQ(ind) = 1.d0;
                Imol(1:Nq) = Imol(1:Nq) + 2.d0*(squeeze(FF(a,b,1:Nq)).').*sinQ(1:Nq); 
            end
        end
        
        if FLAGelec == 0
            Wiam(traj,1:Nq) = Imol(1:Nq) + Iat; % for each traj - add atomic term if X-ray
        else
            Wiam(traj,1:Nq) = (qAng .* Imol(1:Nq)) ./ Iat; %  s*Imol/Iat if electron scattering
        end
    end
    
    %Wiam = 0.25*Wiam; % (1/2me^2)
    Wiam_avg = sum(Wiam,1)./Ntraj; 
    Wiam_tot(1:Nq,ts) = Wiam_avg;
end


if Nts > 1
    pdW = zeros(Nq,Nts);
    for ts=1:Nts
        if FLAGelec == 0 % if not static, calculate percentage diff and plot
            pdW(1:Nq,ts) = exfrac*( Wiam_tot(1:Nq,ts) - Wiam_tot(1:Nq,1) ) ./ Wiam_tot(1:Nq,1);
        else
            pdW(1:Nq,ts) = exfrac*( Wiam_tot(1:Nq,ts) - Wiam_tot(1:Nq,1) );
        end
    end
    [QQ, TT] = meshgrid(qAng, tt);
    mesh(QQ,TT,pdW');
    xlabel(['q (1/' char(197) ')'])
    ylabel('time (fs)')
    xlim([0 7])
else % static plot
    plot(qAng, Wiam_tot);
    xlabel(['q (1/' char(197) ')'])
    xlim([0 7])
end