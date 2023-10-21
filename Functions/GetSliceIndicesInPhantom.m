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
% This obtains the indices of the phantom object that corresponds to the
% current slice. The voxels of the cubes are treated as solid homogenous
% cubes of T1 and T2, therefore, the indices sufficiently inform the T1 and
% T2 values within one slice object.
%
% Inputs:
% SliceIndex = Index of the slice object
% nAngles = Number of angles in the slice windows
% PhantomRes = Resolution of the phantom 
% SliceThickness = Thickness of the slice
% SliceGap = Gap between slices 
% DelZ = Slice Window
% SliceFOV = FOV of the slice window in phantom space
%
% Output:
% SliceIndices = Indices in phantom space of all of the phantom inside the
% slice.

function SliceIndices = GetSliceIndicesInPhantom(SliceIndex,nAngles,PhantomRes,SliceThickness,SliceGap,DelZ,PhantomDims,SliceFOV)

SliceRes = SliceThickness + SliceGap;
PhantomFOV = PhantomDims*PhantomRes;
StartMM = (PhantomFOV - SliceFOV)/2;
SliceCenterPos = ceil(((SliceRes*SliceIndex)+StartMM)/PhantomRes);
SliceIndices = SliceCenterPos+ceil(linspace(-DelZ/PhantomRes,DelZ/PhantomRes,nAngles));

end