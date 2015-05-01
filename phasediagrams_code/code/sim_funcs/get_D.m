function [D,nullvec,D_unmasked] = get_D(mask, bc_type)
%GET_D    Set up finite difference matrix for masked array.
% [D,nullvec,D_unmasked] = get_D(mask, bc_type)
% 
% Inputs:
%  - mask       Square logical array with true meaning pixel is in image.
%  - bc_type    Type of boundary conditions. 'none' means that no boundary
%               conditions are used, and the null space consists of the
%               constant image. 'zero' means that pixels outside the mask
%               are zero, and the null space is trivial.
% Outputs:
%  - D          Finite difference matrix on the pixels in the mask using the
%               specified boundary condition.
%  - nullvec    Matrix where columns form an orthonormal basis for the null
%               space of D.
%  - D_unmasked The finite difference matrix for the square (un-masked)
%               array.
%
% Only works for square masks.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

[s1,s2] = size(mask);
if s1 == s2
    N = s1;
else
    error('First input "mask" must be a square logical array')
end

switch bc_type
    case 'none'
        L = spdiags(ones(N,1)*[-1,1],[0,1],N-1,N);
        I = speye(N);
        D_unmasked = [kron(I,L); kron(L,I)];
        
        % Remove columns for pixels outside mask
        D = D_unmasked(:,mask);
        % Remove rows with only one component left
        D = D(sum(D~=0,2)>1,:);
        
        % Set up null vectors, only if asked for as output
        if nargout > 1
            % For 'none' BC, the constant is the only null vector
            nullvec = ones(size(D,2),1);
            nullvec = nullvec / norm(nullvec);
        end
        
    case 'zero'
        L = spdiags(ones(N,1)*[-1,1],[0,1],N-1,N);
        L(N,N) = -1;
        I = speye(N);
        D_unmasked = [kron(I,L); kron(L,I)];
        
        % Remove columns for pixels outside mask
        D = D_unmasked(:,mask);
        
        % Remove only rows for pixels outside mask
        D = D([mask(:);mask(:)],:);
        
        % Set up null vectors, only if asked for as output
        if nargout > 1
            % For 'zero' BC, the null space is empty
            nullvec = zeros(size(D,2),0);
        end
        
    case 'neumann'
        L = spdiags(ones(N,1)*[-1,1],[0,1],N-1,N);
        L(N,N) = 0;
        I = speye(N);
        D_unmasked = [kron(I,L); kron(L,I)];
        
        % Remove columns for pixels outside mask
        D = D_unmasked(:,mask);
        
        % Remove only rows for pixels outside mask
        D = D([mask(:);mask(:)],:);
        
        % Set to zero the gradient with components outside mask
        D(sum(D~=0,2)==1,:) = 0;
        
        % Set up null vectors, only if asked for as output
        if nargout > 1
            % For 'none' BC, the constant is the only null vector
            nullvec = ones(size(D,2),1);
            nullvec = nullvec / norm(nullvec);
        end
        
    otherwise
        error('unknown bc_type specified')
end
