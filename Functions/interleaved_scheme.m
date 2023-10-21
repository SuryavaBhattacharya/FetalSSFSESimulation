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
% This function builds the Philip's interleaved scheme
% 
% Inputs:
% NbSlices = Number of slices 
%
% Output:
% interleavedSlices_indices = Index of slices following the interleaved
% scheme


function interleavedSlices_indices = interleaved_scheme(NbSlices)

% Input check
if nargin < 1
    error('Missing input(s).');
elseif nargin > 1
    error('Too many inputs.');
end

% Implement interleaved slice acquisition scheme
% if mod(NbSlices,2)==0
% 	interleavedSlices_index = [2:2:NbSlices, 1:2:NbSlices-1];
% else
%     interleavedSlices_index = [2:2:NbSlices-1, 1:2:NbSlices];
% end

interleavedSlices_indices = zeros(1,NbSlices);
StartPoint = 1;
for ii = 1:4
    Nums = ii:4:90;
    EndPoint = length(Nums)+StartPoint-1;
    interleavedSlices_indices(StartPoint:EndPoint) = Nums;
    StartPoint = EndPoint+1;
end


% Display message for debugging
sprintf('Slices will be acquired in the following order:')
fprintf(' %d\n', interleavedSlices_indices(:))

end