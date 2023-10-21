%%% SJM 9-4-15: ISOCHROMAT based Bloch Simulator. Built to have the same syntax
%%% as blochsim_CK which uses CK parameters. This version operates on the magnetization
%%% vectors directly. Can also simulate single vox as long as arrays have 
%%% singleton dimensions correctly defined
%
%   [mxyt,mzt,mxy,mz] = blochsim_isochromat_relax(B1,G,pos,sens,B0,varargin)
%
%   B1  = B1(t,c) = Nt x Nc array
%   G   = [Gx(t) Gy(t) Gz(t)] = Nt x 3 array
%   pos = [x(:) y(:) z(:)] = Ns x 3 array
%   sens= [S1(:) S2(:) ...Sn(:)] = Ns x Nc array
%   B0  = B0(:) = Ns x 1 vector
%
%   Nt = #timepoints, Ns = #space, Nc = #coils
%   B1/B0 are in mT, G in mT/m, pos in m
%
%   args: can include new parameters
%         'dt' followed by new time step in seconds (default 6.4e-6)
%         'M0' followed by new M0 (default [0;0;1])
%         'T1' followed by T1 in seconds (default is inf)
%         'T2' followed by T2 in seconds (default is inf)
%
%   M0 can have dimension 3xNs or 3x1. If 3xNs then pos, sens and B0 must
%   also have the correct dimensions

function [mxy,mz,mxyt,mzt] = blochsim_isochromat_relax(B1,G,pos,sens,B0,varargin)

gam = 267522.1199722082;        % radians per sec per mT
dt = 6.4e-6;
M0=[0;0;1];

T1 = inf;
T2 = inf;

%%% check args
for ii=1:length(varargin)
    if strcmpi(varargin{ii},'dt')
        dt=varargin{ii+1};
    end
    if strcmpi(varargin{ii},'M0')
        M0=varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T1')
        T1=varargin{ii+1};
    end
   if strcmpi(varargin{ii},'T2')
        T2=varargin{ii+1};
    end
 
end

%%% relaxation pars
e1 = exp(-dt/T1);
e2 = exp(-dt/T2);
E = diag([e2 e2 e1]);
E0 = (1-e1); %<-- assumes magnitude of M0 is 1

Ns = size(pos,1);
Nt = size(G,1);


% Temporal evolution of RF field
bxy = sens * B1.';

% Temporal evolution of Bz (due to gradients)
bz = pos * G';

% add off-resonance
bz = bz + repmat(B0(:),1,Nt);

%%% Magnetization 
% Same number of time points as the pulse.
M = zeros([3 Ns Nt]);
% Initialize
if length(M0(:))==3
    M(:,:,1) = repmat(M0,[1 Ns]);%<- start at equilibrium
else
    M(:,:,1) = M0; %<-- user specified a 3xNs array
end
%%% Store 3D rotation matrix in 3x3 matrix, update values rather than
%%% re-define each time (this is slightly faster)
R = eye(3);

%%% Loop over space
for ss = 1:Ns
   
    %%% loop over time
    % Deal with first time point separately
    tt=1;
    rotmat_3d(-gam * dt * [real(bxy(ss,tt));imag(bxy(ss,tt));bz(ss,tt)]);
    M(:,ss,tt) = E*R*M(:,ss,tt);
    M(3,ss,tt) = M(3,ss,tt)+E0;
    for tt = 2:Nt
        rotmat_3d(-gam * dt * [real(bxy(ss,tt));imag(bxy(ss,tt));bz(ss,tt)]);
        M(:,ss,tt,:) = E*R*M(:,ss,tt-1);
        M(3,ss,tt) = M(3,ss,tt)+E0;
    end
end


%%% Outputs
mxyt = squeeze(M(1,:,:) + 1i*M(2,:,:));
mzt = squeeze(M(3,:,:));
    
mxy = mxyt(:,end);
mz = mzt(:,end);

%%% Rotation matrix helper function: declared within scope of parent
%%% function, so R is within scope
function rotmat_3d(u)

% check zero input
if any(u)
    theta=norm(u);
    u=u/theta;
    ct=cos(theta);
    st=sin(theta);
    R(1) = ct + u(1)^2*(1-ct);
    R(2) = u(1)*u(2)*(1-ct)+u(3)*st;
    R(3) = u(1)*u(3)*(1-ct)-u(2)*st;
    R(4) = u(2)*u(1)*(1-ct)-u(3)*st;
    R(5) = ct+u(2)^2*(1-ct);
    R(6) = u(2)*u(3)*(1-ct)+u(1)*st;
    R(7) = u(3)*u(1)*(1-ct)+u(2)*st;
    R(8) = u(2)*u(3)*(1-ct)-u(1)*st;
    R(9) = ct+u(3)^2*(1-ct);
else
    R=eye(3);
end

end

end

