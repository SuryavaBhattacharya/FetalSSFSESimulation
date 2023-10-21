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
% This function checks parameters before running simulation to make sure 
% parameters are possible in a Philip's scanner
% Inputs:
% ACF = Acceleration factor for SENSE
% TE = Echo time
% ESP = Echo spacing
% FOV = Field of View
% ImRes = Resolution of the reconstruction
% Acq = Acquisition matrix
% TSE_Factor = Length of echo train, or number of refocusing pulses
% SliceThickness = Thickness of the slice
% SliceGap = Gap between slices 
% SliceRes = Resolution in Slice direction

function ParamsChecker(ACF,TE,ESP,FOV,ImRes,Acq,TSE_Factor,SliceThickness,SliceGap,SliceRes)

% Check acquisition matrix vs FOV:
Dim = ceil(FOV(1:2)./ImRes);

if Acq < Dim(1)
    msg = "Acquisition matrix too small for given field of view";
    error(msg)
end

% Checks TSE_Factor sufficiently covers k-space:
PastCentre = ACF*(0:(ceil(TE/ESP)-1));
Centre = ceil(Dim(2)/2);
SamplingScheme = [fliplr(Centre+PastCentre),Centre-ACF:-ACF:1];

if length(SamplingScheme) > TSE_Factor
    msg = "Echo train (TSE_Factor) not long enough to cover k-space";
    error(msg)
elseif length(SamplingScheme) < TSE_Factor
    msg = "Echo train (TSE_Factor) too high for given resolution";
    error(msg)
end

% Checks slice thicknes, slice gap and resolution in slice direction
% correspond to each other:
if SliceThickness + SliceGap ~= SliceRes
    msg = "Resolution in slice direction doesn't match slice thickness + " ...
        + "slice gap.";
    error(msg)
end

end