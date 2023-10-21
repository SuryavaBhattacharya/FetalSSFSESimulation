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
% This function samples the k-space for the given slice. It uses the SENSE
% acceleration scheme to do this
%
% Inputs:
% FourierSignals = Time-dependent k-space at each echo.
% ACF = Acceleration factor for SENSE 
% TE = Echo time
% ESP = Echo spacing
% FOV = Field of View 
% ImRes = Resolution of the reconstruction
%
% Output:
% kspace = 2D k-space of current slice

function kspace = Getkspace(FourierSignals,ACF,TE,ESP,FOV,ImRes)
% Find kspace matrix size (this can be different from acquisition matrix
% and while the overall k-space matrix = acquisition matrix, this
% represents the matrix within which data is sampled - the rest of the
% matrix is 0):
Dim = ceil(FOV(1:2)./ImRes);
% Find size of the Fourier representation of each echo:
FOVSignalSize = size(FourierSignals);
% Difference between acquisition matrix and k-space matrix that is sampled
% from:
SizeDiffs = ceil((FOVSignalSize(1:2)-Dim)/2);

% Partial Fourier factor sampling of k-space including SENSE acceleration:
PastCentre = ACF*(0:(ceil(TE/ESP)-1));
Centre = ceil(Dim(2)/2);
SamplingScheme = [fliplr(Centre+PastCentre),Centre-ACF:-ACF:1];
Lines = 1:length(SamplingScheme);

% Sample kspace:
kspace = zeros(FOVSignalSize(1),Dim(2));
for ii = 1:length(SamplingScheme)
    SampLine = SamplingScheme(ii);
    kspace(:,SampLine) = FourierSignals(:,...
        SampLine+SizeDiffs(2),ii);
end

% Fill missing lines with adjacent lines (but, as according to Philip's
% leave the lines outside of the partial Fourier 0):
if ACF ~= 1
    SamplingScheme2 = fliplr(min(SamplingScheme):max(SamplingScheme));
    SamplingScheme2 = setdiff(SamplingScheme2,SamplingScheme);
    for ii = 1:length(SamplingScheme2)
        SampLine = SamplingScheme2(ii);
        Line = max(Lines(SamplingScheme>SampLine));
        kspace(:,SampLine) = FourierSignals(:,...
            SampLine+SizeDiffs(2),Line);
    end
end

% RestOfLines = max(SamplingScheme)+1:Dim(2);
% kspace(:,RestOfLines) = fliplr(conj(kspace(:,FOVSignalSize(2)-RestOfLines+1))')';

% Pad to acquisition matrix size:
kspace = padToSize(kspace,FOVSignalSize(1:2));
end