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
% This function reorients the volume based on Helene Lajous' function in the
% FaBiAN phantom titled "volue_reorient.m"
%
% Inputs:
% Volume = Volume to reorient
% orientation = Orientation (1 - Saggittal, 2 - Coronal, 3 - Axial)
% PutBack = Flag to put the reconstructed stack into the original
% orientation
%
% Output:
% Volume_Reoriented = Reoriented volume

function Volume_Reoriented = volume_reorient2(    Volume, ...
                                             orientation, ...
                                                 PutBack)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs.');
end

Volume_Reoriented = Volume;

% For the sake of clarity, the orientation corresponding to the slice
% thickness will always be given by the 3
switch orientation
    case 1   %sagittal
        if ~PutBack
            Volume_Reoriented = rotate_LabelledImage(Volume,90,[1,0,0],'crop');
        elseif PutBack
            Volume_Reoriented = rotate_LabelledImage(Volume,-90,[1,0,0],'crop');
        end
    case 2   %coronal
        if ~PutBack
            Volume_Reoriented = rotate_LabelledImage(Volume,90,[0,1,0],'crop');
        elseif PutBack
            Volume_Reoriented = rotate_LabelledImage(Volume,-90,[0,1,0],'crop');
        end
end

end