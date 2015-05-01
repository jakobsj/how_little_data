function X = get_test_image(id, N, opts)
%GET_TEST_IMAGE Common interface for creating 2D test images.
% X = get_test_image(id, N, opts)
%
% Inputs:
% - id   : String with name of test image, e.g, 'spikes'.
% - N    : Dimension, generated image is N by N pixels.
% - opts : Struct with options with fields to determine image parameters
%          depending on the type of image as specified by the input id.
%   Possible fields (not all used for all ids):
%   - randstate  : Integer. Random state for rand.
%   - randnstate : Integer. Random state for randn.
%   - k          : Integer. Number of nonzeros in image.
%   - mask       : Logical N-by-N array. Mask specifying which pixels to
%                  consider part of the image, as given by logical 1s in
%                  opts.mask. Pixels corresponding to logical 0s are
%                  considered not part of the produced image.
%   - radius     : Double scalar. Specific to id equal to 'diskspikes' or 
%                  'disksignedspikes'. Specifies radius of small disk 
%                  inside which all nonzeros occur. opts.radius=N/2 
%                  corresponds to the full inscribed disk in the square 
%                  N-by-N array.
%   - intensity  : Double. Pixel value to assign to certain pixel sets for
%                  some ids.
% Output:
% - X : The generated image.
% 
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

switch lower(id)
    case {'spikes', 'signedspikes'}
        rand('state',opts.randstate);
        if strcmpi(id,'spikes')
            vals = rand(opts.k,1);
        else
            vals = 2*rand(opts.k,1) - 1;
        end
        idx  = randperm(sum(opts.mask(:)));
        idx  = sort(idx(1:opts.k));
        x = zeros(sum(opts.mask(:)),1);
        x(idx) = vals;
        X = zeros(N,N);
        X(opts.mask) = x;
    case {'fftpower','fftpower_1_0','fftpower_1_5','fftpower_2_0'}
        randn('state', opts.randnstate);
        d = 2./N;
        x = (-1.+d/2.):d:(1-d/2.);
        [X1,X2] = meshgrid(x);
        amp = (X1.^2 + X2.^2).^(-opts.pp/2.);
        R = randn(N,N);
        phase = R*2.*pi;
        imft = amp.*exp(1j * phase);
        X = abs(ifft2(ifftshift(imft)));
        X = opts.mask .* X;
        [xsort, order]= sort(X(:),'descend');
        xsort(opts.k+1:end) = 0;
        X(order) = xsort;   
    otherwise
        error('Unknown image_id specified')
end
