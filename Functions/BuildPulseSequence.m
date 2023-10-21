% ============================================================================
%
%    Copyright 2023 Suryava Bhattaharya
%    Copyright 2023 King's College London
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%
% ============================================================================
% This function builds the RF flip angle profiles for the excitation and
% train of refocusing pulses.
%
% Inputs:
% exc = Excitation pulse structure
% refoc = Refocusing pulse structure
% TSE_Factor = Length of echo train, or number of refocusing pulses
% delayCollection = Number of preparation pulses
% nAngles = Number of angles in the slice windows
% 
% Output:
% alphas = profile of flip angles across all pulses

function alphas = BuildPulseSequence(exc,refoc,TSE_Factor,delayCollection,nAngles,varargin)
% prepulse = [];
% gam = 2*pi*42577.46778;
% for ii=1:length(varargin)
%     if strcmpi(varargin{ii},'gam')
%         gam = varargin(ii+1);
%     end
%     if strcmpi('prepulse')s
%         prepulse = varargin(ii+1);
%     end
%     if strcmpi('postpulse')
%         postpulse = varargin(ii+1);
%     end
% end

% Variable checker
% outputFigures = false;
% for ii=1:length(varargin)
%     if strcmpi(varargin{ii},'outputFigures')
%         outputFigures = true;
%     end
% end

% Instantiate alphas matrix with all the echoes:
alphas = zeros(nAngles,TSE_Factor+delayCollection);

% Excitation:
alphas(:,1) = (pi/180)*get_rf_angles(exc.rf*exc.amp,nAngles,exc.T,...
    exc.delz,exc.rewinder,exc.G);

% This flag signifies end of preparation pulses:
flag_single = false;

% Going through all the preparation pulses followed by standard refocusing
% pulses:
for ii = 2:TSE_Factor+delayCollection+1
    if flag_single
        alphas(:,ii) = alphas(:,ii-1);
    end
    
    % Assign profile from angles: 
    alphas(:,ii) = (pi/180)*get_rf_angles(refoc.amp(1)*refoc.rf,nAngles,...
        refoc.T(1),refoc.delz(1),...
        refoc.rewinder(1),refoc.G(1));
    
    % Deal with the amplitude and duration arrays (remove used amplitudes
    % or durations until only one is left for the rest of the train):
    if length(refoc.amp)==1
        outputFigures = false;
    end
    if length(refoc.amp)>1
        refoc.amp = refoc.amp(2:end);
    end
    if length(refoc.T)>1
        refoc.T = refoc.T(2:end);
    end
    if length(refoc.delz)>1
        refoc.delz = refoc.delz(2:end);
    end
    if length(refoc.rewinder(1))>1
        refoc.rewinder = refoc.rewinder(2:end);
    end
    if length(refoc.G(1))>1
        refoc.G = refoc.G(2:end);
    end

    % When a single refocusing pulse amplitude and duration left, then
    % insure that no more amplitudes and durations are deleted:
    if ~flag_single && length(refoc.amp) == 1 && length(refoc.T) == 1 ...
            && length(refoc.delz) == 1 && length(refoc.rewinder(1)) == 1 ...
            && length(refoc.G(1)) == 1
        flag_single = true;
    end
end
end

function alphas = get_rf_angles(pulseRF,nAngles,T,delz,rewinder,G,varargin)
% Helper function to get flip angle profiles

% Get magnetisation profile:
[~,~,~,mzt]= get_sliceProfile(pulseRF,T,delz,rewinder,G,nAngles);

% Get magnetisation at start
M_vec_t_1 = [zeros(size(mzt(:,1)'));mzt(:,1)'];

% Get magnetisation at the end of the pulse
M_vec_t_end = [sqrt(mzt(:,1).^2-mzt(:,end).^2)';mzt(:,end)'];

% Get flip angles from start magnetisation and end magnetisation.
CosTheta = dot(M_vec_t_end,M_vec_t_1)./(vecnorm(M_vec_t_1).*vecnorm(M_vec_t_end));
alphas = acosd(CosTheta);
% degProf = acosd(CosTheta);
% degProf = 90-rad2deg(atan(mz./real(mxy)));
% midpnt = ceil(length(degProf)/2);
% firstpart = linspace(1,midpnt,ceil(nAngles/2));
% secondpart = linspace(midpnt+firstpart(2)-firstpart(1),length(degProf),floor(nAngles/2));
% inds = [firstpart,secondpart];
% alphas = degProf(round(inds));

end

function [mxy,mz,mxyt,mzt]= get_sliceProfile(rf,T,delz,rewinder,Gsel,nAngles)
% Helper function to prepare pulses for the Bloch Simulations
% Assign time intervals:
dt = T/length(rf);
% Assign number of time points:
nt = length(rf);
% Deal with Gradients:
G = zeros([nt 3]);
G(:,3)=Gsel;
% Deal with rewinder:
if rewinder==true
    rewinder_G = zeros([floor(nt/2) 3]);
    rewinder_G(:,3)= -1*Gsel;
    G = [G; rewinder_G];
    rf_rewinder = zeros([floor(nt/2), 1]);
    rf = [rf, rf_rewinder'];
%     nt = ceil(1.5*nt);
%     nAngles = ceil(1.5*nAngles);
end
% Spatial points:
zz=linspace(-1*delz,delz,nAngles);
% Assign the magnetisation vectors across the slice
pos = zeros([nAngles 3]);
pos(:,3)=zz(:);
%%% Simulate using Bloch simulations
[mxy,mz,mxyt,mzt] = blochsim_isochromat(rf',G,pos,ones(nAngles,1),zeros(nAngles,1),'dt',dt);
end