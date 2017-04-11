%KTRUSS Computer the k-truss of a graph
% Inputs:
%       inc_mat_file - path to file containing incidence matrix
%       k 
%
function [E,t_ktruss] = ktruss(inc_mat_file, k)

ii = load(inc_mat_file);
t0 = clock;
E = sparse( ii(:,1), ii(:,2), ii(:,3) );
t_read = etime(clock, t0);
fprintf(2, 'time to read data  : %f sec.\n', t_read);

fprintf(2, 'computing k-truss\n');
t0 = clock;
tmp = E.'*E;
R = E* (tmp-diag(diag(tmp)));

s = sum(double(R==2),2);
xc = s >= k-2;
[m,n] = size(E);
% While edges exist violating k-Truss, delete those edges and take a subgraph.
while nnz(xc) ~= nnz(any(E,2)) 
    Etrans = E.';
    Ekeep = Etrans(:, xc);
    [i,j,v] = find(Ekeep);
    E = sparse(find(xc)(j), i, v, m, n);
    %E(not(xc),:) = 0;
    tmp = E.'*E;
    R = E* (tmp-diag(diag(tmp)));
    s = sum(double(R==2),2);
    xc = s >= k-2;
end

t_ktruss = etime(clock, t0);
fprintf(2, 'time to compute k-truss : %f\n', t_ktruss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph Challenge benchmark
% Developer : Dr. Jeremy Kepner (kepner@ll.mit.edu)
%             Dr. Siddharth Samsi (ssamsi@mit.edu)
%
% MIT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) <2017> Massachusetts Institute of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

