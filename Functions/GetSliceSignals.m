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
% SliceObject = cell of 3D T1 and T2map of the spatial distribution of T1
% and T2 in the map.
% Alphas = Flip angles along the slice
% Echo spacing
%
% Output:
% SliceSignals = Spatial profile of the current slice

function SliceSignals = GetSliceSignals(SliceObject,Alphas,ESP)
%tic

% Find size of Slice array
Dims = size(SliceObject{1});
Dims = Dims(1:2);
% Instantiate Signal array:
SliceSignals = zeros([size(Alphas,1),Dims,size(Alphas,2)-1]);

% To correct for slice profile, loop through all the positions on the slice
% window:
for aa = 1:size(Alphas,1)
    % Get the 2D T1 and T2 map of the current position in the slice window
    T1map = squeeze(SliceObject{1}(:,:,aa));
    T2map = squeeze(SliceObject{2}(:,:,aa));
    % Skip through 0 values
    if sum(T1map(:)) == 0
        continue;
    end
    % Flatten the map so it's easier to deal with:
    [Vals, inds] = unique([T1map(:), T2map(:)], 'rows');
    SliceSecSignals = zeros([length(T1map(:)),size(Alphas,2)-1]);
    % Go through all the pixels: 
    for ii = 1:length(inds)
        T1 = Vals(ii,1);
        T2 = Vals(ii,2);
        if T1 == 0 
            continue
        end
        % EPG Simulation:
        SliceSecSignals(T1map(:)==T1,:) = repelem(EPG_TSE(Alphas(aa,:),ESP,...
            T1,T2,'kmax',size(Alphas,2)+2),length(T1map(T1map==T1)),1);
        
    end
    % Reorganise into the shape of the map
    SliceSignals(aa,:,:,:) = reshape(SliceSecSignals,[Dims,size(Alphas,2)-1]);
end
% Sum/Integrate over slice window.
SliceSignals = squeeze(sum(SliceSignals,1));
%tSignalGen = toc;
%disp("Time taken to generate EPG signals: " + num2str(tSignalGen) + "s")
end
