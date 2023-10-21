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
% This function obtains the T1 and T2 values in the slice, dealing with the
% slice as an effective 3D volume spanning the width of the slice window.
%
% Inputs:
% T1Phantom = T1 map of the phantom
% T2Phantom = T2 map of the phantom
% SliceIndices = Indices in phantom space of all of the phantom inside the
% slice.
%
% Output:
% SliceObject = cell of 3D T1 and T2map of the spatial distribution of T1
% and T2 in the map.


function SliceObject = GetSliceObject(T1Phantom,T2Phantom,SliceIndices)

SliceObject = cell(1,2);
SliceObject{1} = T1Phantom(:,:,SliceIndices);
SliceObject{2} = T2Phantom(:,:,SliceIndices);


end