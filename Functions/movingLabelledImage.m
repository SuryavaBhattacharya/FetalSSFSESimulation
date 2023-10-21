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
% This function applies transformations to the labelled segmentation image 
% and uses linear interpolation on individual tissue classes followed by 
% max voting to perform better transformations
%
% Inputs:
% LabelledImage = Map of tissue segmentations (phantom)
% ImRes = Resolution of the reconstruction
% translation_displacements = displacements for translation in mm in 3D
% rotation_angles = Angles for rotation in ° in any direction in 3D 
% rotation_axes = Axes among which rotations were made
% motion_level = Level of motion
% bbox = State of matrix after transformation
%
% Output:
% movedLabelledImage = Transformed segmentation maps
% rotationDone = Logical informing if a rotation was actually made

function [movedLabelledImage,rotationDone] = movingLabelledImage(LabelledImage,...
                                                                ImRes, ...
                                             translation_displacement, ...
                                                       rotation_angle, ...
                                                        rotation_axis, ...
                                                         motion_level, ...
                                                                 bbox)
% Find all tissue labels:
uniqueLabels = unique(LabelledImage);
rotationDone = false;
% Instantiate random sampler: 
sampler = true(1,2);
switch motion_level
    case 0  %no motion
    case 1  %slight motion 10% motion chance
        sampler = [true,false(1,9)];
    case 2  %moderate motion 12.5% motion chance
        sampler = [true,false(1,7)];
    case 3  %strong motion 16.7% motion chance
        sampler = [true,false(1,5)];
    case 4  %extreme motion 25% motion chance
        sampler = [true,false(1,3)];
end

% Input check
if nargin < 7
    error('Missing input(s).');
elseif nargin > 7
    error('Too many inputs.');
end

% Check random sampling chance
if randsample(sampler,1)
    tic

    Reference_LabelledImage_mm = imref3d(size(LabelledImage), ImRes, ImRes, ImRes);
    % Keep a reference image for max voting:
    movedLabelledImage = zeros(size(LabelledImage));

    % Probability map for max voting:
    LastLabelProbMap = zeros(size(LabelledImage));

    % Loop through tissue classes
    for ii = 1:length(uniqueLabels)
        % Mask current segmentation:
        BinMaskLab = double(LabelledImage==uniqueLabels(ii));
        % Translate:
        BinMaskLab_translated= imtranslate(BinMaskLab, ...
            Reference_LabelledImage_mm, translation_displacement,  'linear');

        % Apply 3D rotation, defined by a rotation angle and a rotation axis, to
        % the already translated fetal brain anatomy:
        BinMaskLab_rotated = imrotate3(BinMaskLab_translated, rotation_angle, rotation_axis, 'linear', bbox);

        % Update segmentation map based on maximum probability (max voting):
        movedLabelledImage(BinMaskLab_rotated>LastLabelProbMap) = uniqueLabels(ii);
        % Update probability map:
        LastLabelProbMap = reshape(max([BinMaskLab_rotated(:),LastLabelProbMap(:)],[],2),size(movedLabelledImage));
    end
    tRotate = toc;
    disp("Time taken to rotate labelled image " + num2str(tRotate) + "s.")
    rotationDone = true;
else
    % If no motion to be applied keep original tissue:
    movedLabelledImage = LabelledImage;
end

end