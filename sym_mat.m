function mat = sym_mat(name, m, n)

if(nargin == 3)
    dim = [m n];
elseif (nargin == 2)
    dim = m;
else
    error('wrong number of arguments!')
end
        

mat = sym(name, dim);
assume(mat, 'real');