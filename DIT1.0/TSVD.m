function [U,S,V]=TSVD(mat,tol)

    %  This fuction calculates Truncated Singular Value Decomposition (TSVD)
    %  [U,S,V]=TSVD(matrix,Tolerance)
    
    if nargin==1
        tol=eps;
    end
    
    [Ut,St,Vt]=svd(mat);
    Sdia=diag(St);
    box=find((Sdia>=tol),1,'last');
    [U,S,V]=svds(mat);

end