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
% This function applies rotations to the labelled segmentation image
% and uses linear interpolation on individual tissue classes followed by
% max voting to perform better transformations
%
% Inputs:
% LabelledImage = Map of tissue segmentations (phantom)
% rotation_angles = Angles for rotation in ° in any direction in 3D 
% rotation_axes = Axes among which rotations were made
% motion_level = Level of motion
% bbox = State of matrix after transformation
%
% Output:
% rotatedLabelledImage = Transformed segmentation maps

function rotatedLabelledImage = rotate_LabelledImage(LabelledImage,rotation_angle,rotation_axis,bbox)

uniqueLabels = unique(LabelledImage);

rotatedLabelledImage = zeros(size(LabelledImage));

LastLabelProbMap = zeros(size(LabelledImage));

for ii = 1:length(uniqueLabels)
    BinMaskLab = double(LabelledImage==uniqueLabels(ii));

    % Apply 3D rotation, defined by a rotation angle and a rotation axis, to
    % the already translated fetal brain anatomy
    BinMaskLab_rotated = imrotate3(BinMaskLab, rotation_angle, rotation_axis, 'linear', bbox);

    rotatedLabelledImage(BinMaskLab_rotated>LastLabelProbMap) = uniqueLabels(ii);
    LastLabelProbMap = reshape(max([BinMaskLab_rotated(:),LastLabelProbMap(:)],[],2),size(rotatedLabelledImage));
end
end