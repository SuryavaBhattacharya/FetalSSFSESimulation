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
% This function creates a k-space representation of the slice of the 
% spatial distribution of signal intensities at each time point over the 
% slice window. This deals with zero padding and cropping of k-space to 
% deal with the spatial resolution and FOV.
%
% Inputs:
% SliceSignals = Spatial profile of the current slice
% FOV = Field of View 
% PhantomRes = Resolution of the phantom 
% Acq = Acquisition matrix
%
% Output:
% FourierSignals = Time-dependent k-space at each echo.


function FourierSignals = ConvertSignalsToFourier(SliceSignals,FOV,PhantomRes,Acq)
%tic
% We first obtain the appropriate FOV in spatial domain and crop in (or 
% zero pad out). Then we deal with resolution and acquisition matrix in
% k-space, which simulates how the scanner deals with k-space. 
% NOTE: This only works when using a high resolution phantom and down
% sampling. In order to upsample, one would have to upsample the phantom
% itself, as a first step, to the higher resolution before any simulating

% Get FOV in slice:
FOV = FOV(1:2);

% Get the FOV covered by the phantom matrix size:
SignalsDims = size(SliceSignals);
PhantomFOV = PhantomRes(1:2).*SignalsDims(1:2);
% The target we want to crop/pad the signals array to:
TargetDims = ceil(FOV./PhantomRes(1:2));
% Pad or crop:
if PhantomFOV<FOV
    PadSize = ceil((TargetDims-SignalsDims(1:2))/2);
    SignalsCopy = padarray(SliceSignals,[PadSize,0]);
    if size(SignalsCopy,1) ~= TargetDims(1)
            SignalsCopy = SignalsCopy(1:end-1,:,:);
    end
    if size(SignalsCopy,2) ~= TargetDims(2)
            SignalsCopy = SignalsCopy(:,1:end-1,:);
     end
    SliceSignals = SignalsCopy;
    clear SignalsCopy;
elseif PhantomFOV > FOV
    StartPosition = floor((SignalsDims(1:2)-TargetDims)/2)+1;
    SliceSignals = SliceSignals(StartPosition(1):StartPosition(1)+TargetDims(1)-1,...
        StartPosition(2):StartPosition(2)+TargetDims(2)-1,:);
end

% The Fourier part:
FourierSignals = zeros([TargetDims,SignalsDims(3)]);

% Loop through all echoes and FFT the spatial distribution of signal
% intensities:
for ii = 1:SignalsDims(3)
    FourierSignals(:,:,ii) = ...
        fftshift(fft2(squeeze(SliceSignals(:,:,ii))));
end

% Crop to the appropiate acuisition matrix:
StartPosition = floor((TargetDims-Acq)/2)+1;
FourierSignals = FourierSignals(StartPosition(1):StartPosition(1)+Acq-1,...
    StartPosition(2):StartPosition(2)+Acq-1,:);

%tConvertFourier = toc;
%disp("Time taken to convert to Fourier: " + num2str(tConvertFourier) + "s")
end