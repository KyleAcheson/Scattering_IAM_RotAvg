% Kyle Acheson - 07/09/2010
% IAM rotationally averaged 1D scattering pattern calculator
% NOW loops over Ntraj & Nts - NOW sum trajectories for tot. I(q) 
% Uses Adams f_functions code for form factors - H,C,N,F,S,I,Xe available.
%clear all 
%close all 

tstart = tic;

%load('sf6xyz.mat'); % Reads distance matrix (Natom,Natom,Ntrj,Nts)

%Q = geom; % must be in angstroms
[x,y,Ntraj,Nts] = size(Q);

% Setup
Ang=char(197);
FLAGelec = 0; % 0 for x-ray, 1 for UED 
au2ang = 0.52917721092d0;
ang2au = 1/au2ang;
Q = Q%*ang2au;
atmnum = [54]; % atoms
kin   = 12.d0*au2ang; % incident wave vector inv au
theta_lim = [0 pi]; % range of theta angles
Nq = 481; % number of points over theta range 
Nphi = 1; % number of phi - only theta matters in rot. average 

theta_min = theta_lim(1);
theta_max = theta_lim(2);
qmin =2.d0*kin*sin(0.5d0*theta_min);
qmax = 2.d0*kin*sin(0.5d0*theta_max);
q = linspace(qmin,qmax,Nq); % linear in q
qAng = q*ang2au; % q in inv Ang



[FF,fq]  = get_scattering_factors(qAng,atmnum,FLAGelec); % NB: call using q in Angstroms^(-1)
Iat      = sum(fq.^2); % atomic scattering term 

Wiam_tot = zeros(Ntraj,Nq,Nts);
for ts=1:Nts % loop over trajs and time steps
    Wiam = zeros(Ntraj,Nq);
    Q = Q(:,:,:,ts);
    for traj=1:Ntraj
        tot = zeros(1,Nq);
        sinQ = zeros(1,Nq);
        for a=1:x
            for b=a+1:x
                DQ = q(1:Nq)*norm(Q(a,1:3,traj)-Q(b,1:3,traj)); %qRij
                sinQ(1:Nq) = sin(DQ(1:Nq))./DQ(1:Nq); % Sin(qRij)/qRij
                ind = find(abs(DQ)<1.d-9); % check q=0 limit
                sinQ(ind) = 1.d0;
                tot(1:Nq) = tot(1:Nq) + 2.d0*(squeeze(FF(a,b,1:Nq)).').*sinQ(1:Nq);
            end
        end
        Wiam(traj,1:Nq) = tot(1:Nq) +Iat; % add atomic term
    end
    Wiam_tot(1:Ntraj,1:Nq,ts) = Wiam; % total scattering for all trajs
end
Wiam_tot = Wiam_tot*0.25; % (e/2m_e)^2 in a.u.
Wiam_avg = sum(Wiam_tot,1)./Ntraj; % Ehrenfest wf scattering - equal weights = avg.