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
% This function pads any image to a target size
%
% Inputs:
% Image = Image to pad
% targetSize = Size of target image (up to 4D image)
%
% Output:
% paddedImage = Image padded to the right size

function paddedImage = padToSize(Image,targetSize)

% Get image dimension and pad size. For an odd difference we might need to
% remove one of the rows/colums etc to get to the right size:
Dims = size(Image);
PadSize = zeros(1,length(Dims));
ToRemove = false(1,length(Dims));

% Get pad size and check if the difference between target size in dimension
% and image size is even or odd to remove an element in said dimension:
for ii = 1:length(Dims)
    Diffam = targetSize(ii) - Dims(ii);
    if mod(Diffam,2)~= 0
        ToRemove(ii) = true;
    end
    PadSize(ii) = ceil(Diffam/2);
end

% Pad array:
paddedImage = padarray(Image,PadSize);

% Remove elements if necessary in appropriate dimensions:
if ToRemove(1)
    paddedImage = paddedImage(1:end-1,:,:,:);
end
if ToRemove(2)
    paddedImage = paddedImage(:,1:end-1,:,:);
end
if length(ToRemove) ==3
    if ToRemove(3)
        paddedImage = paddedImage(:,:,1:end-1,:);
    end
end

if length(ToRemove) ==4
    if ToRemove(4)
        paddedImage = paddedImage(:,:,:,1:end-1);
    end
end

end