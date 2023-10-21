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
% This function applies scaling to the labelled segmentation image
% and uses linear interpolation on individual tissue classes followed by
% max voting to perform better transformations
%
% Inputs:
% LabelledImage = Map of tissue segmentations (phantom)
% factor = Factor to scale image by
%
% Output:
% rotatedLabelledImage = Transformed segmentation maps

function ResampledLabelledImage = ResizedLabelledImage(LabelledImage,factor)

uniqueLabels = unique(LabelledImage);

ResampledLabelledImage = zeros(size(LabelledImage));



for ii = 1:length(uniqueLabels)
    BinMaskLab = double(LabelledImage==uniqueLabels(ii));

    % Apply 3D rotation, defined by a rotation angle and a rotation axis, to
    % the already translated fetal brain anatomy
    if (isequal(round(factor),factor) && isempty(factor(factor<10))) || length(factor) == 1
        BinMaskLab_resized = imresize3(BinMaskLab,factor, 'Method', 'linear');
    else
        BinMaskLab_resized = imresize3(BinMaskLab,'Scale', factor, 'Method', 'linear');
    end
    if ii == 1
        LastLabelProbMap = BinMaskLab_resized;
        ResampledLabelledImage = zeros(size(BinMaskLab_resized));
    else
        ResampledLabelledImage(BinMaskLab_resized>LastLabelProbMap) = uniqueLabels(ii);
        LastLabelProbMap = reshape(max([BinMaskLab_resized(:),LastLabelProbMap(:)],[],2),size(ResampledLabelledImage));
    end
end


end