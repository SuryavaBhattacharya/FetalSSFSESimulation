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
% This function is based on Helene Lajous' "motion_transform.m" in the
% FaBiAN phantom to obtain motion parameters by random from a normal 
% distribution.
%
% Inputs:
% motion_level = Level of motion
% NbSlices = Number of slices 
% PhantomRes = Resolution of the phantom 
%
% Output:
% translation_displacements = displacements for translation in mm in 3D
% rotation_angles = Angles for rotation in ° in any direction in 3D 
% rotation_axes = Axes among which rotations were made

function [translation_displacements, ...
                    rotation_angles, ...
                      rotation_axes] = GetMotionTransforms(motion_level, ...
                                                               NbSlices, ...
                                                             PhantomRes)

% Input check
if nargin < 3
    error('Missing input(s).');
elseif nargin > 3
    error('Too many inputs.');
end

% We define 3 levels of motion where 5% of the total number of slices is
% affected by 3D random rigid motion of the fetus:
%  1 - Slight motion:   translation will vary between [-1,+1]mm and
%                       rotation between [-2,+2]°;
%  2 - Moderate motion: translation will vary between [-3,+3]mm and
%                       rotation between [-5,+5]°;
%  3 - Severe motion:   translation will vary between [-6,+6]mm and
%                       rotation between [-10,+10]°.
%  4 - Extreme motion:   translation will vary between [-12,+12]mm and
%                       rotation between [-20,+20]°.
switch motion_level
    case 0  %no motion
    case 1  %slight motion
        translation_amplitude = 1;
        rotation_amplitude = 2;
    case 2  %moderate motion
        translation_amplitude = 3;
        rotation_amplitude = 5;
    case 3  %strong motion
        translation_amplitude = 6;
        rotation_amplitude = 10;
    case 4  %extreme motion
        translation_amplitude = 12;
        rotation_amplitude = 20;
end

translation_amplitude = translation_amplitude*PhantomRes;
% Initialization of tranformation arrays
translation_displacements = zeros(NbSlices, 3);
rotation_axes = zeros(NbSlices, 3);
rotation_angles = zeros(NbSlices, 1);

% Loop through slices
for iSlice=1:NbSlices
    % All parameters are normally distributed: 
    % Translation
    translation_displacements(iSlice,:) = normrnd(0,(1/3)*translation_amplitude,1,3); %in mm
    % Rotation
    rotation_angles(iSlice) = normrnd(0,(1/3)*rotation_amplitude); %in degrees
    rotation_axes(iSlice,:) = [unifrnd(0,1) unifrnd(0,1) unifrnd(0,1)]; 
    % Display message for debugging
    sprintf('The motion corrupted slice #: %d is characterized by a translational displacement of [%0.4f %0.4f %0.4f]mm along the 3 main axes, and by a rotation of %0.4f degrees around the axis [%0.4f %0.4f %0.4f].', iSlice, translation_displacements(iSlice,:), rotation_angles(iSlice), rotation_axes(iSlice,:))
end

end