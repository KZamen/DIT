function [nlow,nhigh]=counter(S,tol)
%%% This function finds number of diagonal elemnts of a matrix that are
%%% larger than a given tolerance.
%%% [nlow,nhigh]=counter(S,tol)
%%% nlow: Number of diagonal elements that are smaller than tolerance
%%% nhigh: Number of diagonal elements that are bigger than tolerance
%%% S: investigated matrix
%%% tol: Tolerance

    keep=diag((S)>tol);
    keep=sort(keep,'descend');
    numb=find(keep>0);
    nhigh=size(numb,1);
    nlow=size(diag(S),1)-nhigh;
end