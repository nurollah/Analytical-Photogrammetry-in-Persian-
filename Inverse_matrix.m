% Bismillah
function [A_inv]=Inverse_matrix(A)
[m,n,h]=size(A);
A_inv = zeros(m,n);
if m~=n
    disp("this matrix could not be inveresed");
else
    
    det1 = det(A);
    if det1==0
        disp("the determinat of input matrix is zero.");
    else
        for i=1:m
            for j=1:n
                Mij = A;
                Mij(i,:)=[];
                Mij(:,j)=[];
                Cji=(-1)^(i+j)*det(Mij);
                A_inv(j,i)= Cji / det1;
            end
        end
    end
end
end

