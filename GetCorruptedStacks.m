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
% This function builds the motion corrupted stacks from a single-shot fast
% spin-echo sequence based on a Philip's system
%
% Inputs:
% Fetal_Brain_model_path = Path to phantom
% GA = Gestational age 
% exc = Excitation pulse structure
% refoc = Refocusing pulse structure
% nAngles = Number of angles in the slice windows
% NbSlices = Number of slices 
% SliceWindow = The window to sample signals related to the slice
% ImRes = Resolution of the reconstruction
% Acq = Acquisition matrix
% TSE_Factor = Length of echo train, or number of refocusing pulses
% FOV = Field of View 
% B0 = Field strength of system
% TE = Echo time
% ESP = Echo spacing
% ACF = Acceleration factor for SENSE 
% motion_level = Level of motion
% PhantomRes = Resolution of the phantom 
% SliceThickness = Thickness of the slice
% SliceGap = Gap between slices 
% std_noise = standard deviation of Gaussian noise distribution
% orientation = Orientation (1 - Saggittal, 2 - Coronal, 3 - Axial)
% StackName = Output name of the motion corrupted stacks
%
% Output:
% CorruptedStack = Stack of motion corrupted stacks - also written to a
% nifti file

function CorruptedStack = GetCorruptedStacks(Fetal_Brain_model_path,GA,exc,refoc,nAngles,NbSlices, ...
    SliceWindow,ImRes,Acq,TSE_Factor,FOV,B0,TE,ESP,ACF,motion_level,PhantomRes, ...
    SliceThickness,SliceGap,std_noise,orientation,StackName)

% This function generates a stack of motion corrupted slices with
% a single-shot fast spin echo sequence

% Add subpaths:
addpath("bin/","EPGX-src/","lib/","Functions/")

% Verify parameters and throw errors if improper parameters used:
ParamsChecker(ACF,TE,ESP,FOV,ImRes,Acq,TSE_Factor,SliceThickness,SliceGap,ImRes);

% Can simulate no motion:
SimulateMotion = true;
if motion_level < 1
    SimulateMotion = false;
end

% Obtain flip angle profile over the sequence:
Alphas =  BuildPulseSequence(exc,refoc,TSE_Factor,0,nAngles);

% Read in phantom
Fetal_Brain = niftiread(strcat(...
    Fetal_Brain_model_path, 'STA', sprintf('%02s', ...
    num2str(GA)), '_tissue.nii.gz'));

% Pad array for transformation to avoid losing parts of the brain in
% transformations:
Fetal_Brain = padarray(Fetal_Brain,[10,10,10]);

% Reorient to desired view for slice sampling:
Fetal_Brain = volume_reorient2(Fetal_Brain,orientation,false);
reorientedDims = size(Fetal_Brain);

% Consider only non-zero labels:
Fetal_Brain_Labels = unique(Fetal_Brain);
Fetal_Brain_Labels = Fetal_Brain_Labels(Fetal_Brain_Labels > 0);
% Write labels of the different segmented brain tissues in a list:
Fetal_Brain_Tissues = permute(Fetal_Brain_Labels, [2,1]);


% Philip's interleaved scheme:
interleavedSlices_indices = interleaved_scheme(NbSlices);

delz = SliceWindow;

%Identify rotated slices:
rotatedSlices = [];

% Obtain transformed brain for all of the sampled slices and avoid storing
% a large signal array in memory:
RotatedPhantom = zeros([NbSlices,size(Fetal_Brain)]);
if SimulateMotion
    [translation_displacements, ...
        rotation_angles, ...
        rotation_axes] = GetMotionTransforms(motion_level, ...
        NbSlices,PhantomRes(1));

    % For light and moderate motion, the matrix can be cropped, to speed
    % up the transformations
    if ~ismember(motion_level,[3,4])

        % if the brain is too big, the matrix can't be cropped this try 
        % catch statement avoids cropping into brains that are too big.
        try
            [rows,columns,slices] = ind2sub(size(Fetal_Brain),find(Fetal_Brain));

            Fetal_Brain_2_rotate = Fetal_Brain(min(rows)-6:max(rows)+6,...
                min(columns)-6:max(columns)+6,min(slices)-6:max(slices)+6);
        catch
            % Brains too big, just use uncropped matrix:
            Fetal_Brain_2_rotate = Fetal_Brain;
        end
    else
        % High motion brains, just use uncropped matrix
        Fetal_Brain_2_rotate = Fetal_Brain;
    end

    % Instantiate matrix of transformed brains
    RotatedPhantom = zeros([NbSlices,size(Fetal_Brain_2_rotate)]);
    
    % Store transformations for all slice instances.
    for ii = 1:NbSlices
        [Fetal_Brain_2_rotate, rotationDone] = movingLabelledImage(Fetal_Brain_2_rotate, ...
            PhantomRes(1), ...
            translation_displacements(ii,:), ...
            rotation_angles(ii), ...
            rotation_axes(ii,:), ...
            motion_level, ...
            'crop'); % < Store previous transformation as the fetus doesn't move back
        RotatedPhantom(ii,:,:,:) = Fetal_Brain_2_rotate;
        if rotationDone
            rotatedSlices = [rotatedSlices, interleavedSlices_indices(ii)];
        end
    end

end

% Make sure the size of the phantom matric is sufficient to sample the
% slices with the relevant slice windows:
if ceil(PhantomRes(3)*reorientedDims(3)) < ceil(NbSlices*ImRes+2*SliceWindow) + 2
    Fetal_Brain = padToSize(Fetal_Brain,[reorientedDims(1:2),ceil(((NbSlices*ImRes+2*SliceWindow) + 2)/PhantomRes(3))]);

    if SimulateMotion
        TargetRotatedPhantomSize = [size(RotatedPhantom,1),reorientedDims(1:2),...
            ceil(((NbSlices*ImRes+2*SliceWindow) + 2)/PhantomRes(3))];
        RotatedPhantom = padToSize(RotatedPhantom,TargetRotatedPhantomSize);
    end
else
    if SimulateMotion
        TargetRotatedPhantomSize = [size(RotatedPhantom,1),reorientedDims];
        RotatedPhantom = padToSize(RotatedPhantom,TargetRotatedPhantomSize);
    end
end

% Generate reference T1 and T2 maps of the fetal brain (for no motion
% case):
[T1map, T2map] = tissue_to_MR(Fetal_Brain, Fetal_Brain_Tissues, B0);

% Loop through slices in parallel loop:
Stacks = zeros(Acq,Acq,NbSlices);
parfor ii = 1:NbSlices
    % Instantiate maps:
    T1mapMoving = T1map;
    T2mapMoving = T2map;
    Slice_ii= interleavedSlices_indices(ii);
    if SimulateMotion
        % Generate reference T1 and T2 maps of the fetal brain:
        [T1mapMoving, T2mapMoving] = tissue_to_MR(squeeze(RotatedPhantom(ii,:,:,:)), Fetal_Brain_Tissues, B0);
    end

    % Build 3D slice object:
    SliceIndices = GetSliceIndicesInPhantom(Slice_ii,nAngles,PhantomRes(3),SliceThickness,SliceGap,...
        delz,size(T1mapMoving,3),FOV(3));
    SliceObject = GetSliceObject(T1mapMoving,T2mapMoving,SliceIndices);

    % Obtain signals and their k-space representations: 
    SliceSignals = GetSliceSignals(SliceObject,Alphas,ESP);
    FourierSignals = ConvertSignalsToFourier(SliceSignals,FOV,PhantomRes,Acq);
    KSpace = Getkspace(FourierSignals,ACF,TE,ESP,FOV,ImRes);
    KSpaceNoise = add_noise(KSpace,std_noise);

    % Obtain spatial version of the slice:
    Stacks(:,:,ii) = ifft2(KSpaceNoise);
end

% Fix slice ordering: 
StacksCopy = zeros(Acq,Acq,NbSlices);
for ii = 1:NbSlices
    Slice_ii= interleavedSlices_indices(ii);
    StacksCopy(:,:,Slice_ii) = Stacks(:,:,ii);
end

% Convert to int16 for Nifti output: 
%CorruptedStack = int16(100*abs(volume_reorient2(StacksCopy,orientation,true)));
CorruptedStack = int16(100*abs(StacksCopy));

View = ["Sagittal", "Coronal", "Axial"];

% Rename stack with view and write image:
StackName = StackName + View(orientation);

niftiwrite(CorruptedStack,...
    StackName    , ...
    'Compressed',true)

% Fix header and orientation for SVRTK and IRTK:
Infothing = niftiinfo(StackName);

Infothing.TransformName = 'Qform';
Infothing.ImageSize = size(CorruptedStack);
Infothing.Datatype = 'int16';
Infothing.raw.pixdim = [1,Infothing.PixelDimensions,1,0,0,0];
switch orientation
    case 1
        Infothing.PixelDimensions = [ImRes,FOV(1:2)/Acq];
        Infothing.Transform.T(1:3,1:3) = Infothing.Transform.T(1:3,1:3)*roty(90)*rotz(180);
    case 2
        Infothing.PixelDimensions = [FOV(1)/Acq,ImRes,FOV(2)/Acq];
        Infothing.Transform.T(1:3,1:3) = Infothing.Transform.T(1:3,1:3)*rotx(90)*rotz(180);
    case 3
        Infothing.PixelDimensions = [FOV(1:2)/Acq,ImRes];
        Infothing.Transform.T(1:3,1:3) = Infothing.Transform.T(1:3,1:3)*rotz(180);
end

% Rewrite image with correct orientation for SVRTK and IRTK:
niftiwrite(CorruptedStack,...
    StackName, Infothing, ...
    'Compressed',true)

% Only works if MIRTK or SVRTK installed (puts images on 0 origin):
cmd = "edit-image " + StackName + ".nii.gz " +...
    StackName + ".nii.gz " + ...
    "-origin 0 0 0";
system(cmd)

end
