
inc_mtx_file = '../../../data/amazon0302_inc.tsv';

E_expected =  [1  1  0  0  0; ...
               0  1  1  0  0; ...
               1  0  0  1  0; ...
               0  0  1  1  0; ...
               1  0  1  0  0; ...
               0  0  0  0  0];


E = ktruss(inc_mtx_file, 3);

[rowinds, colinds, nzvals] = find(E);
dlmwrite('amazon0302_result.tsv', [rowinds colinds nzvals], '\t');

%if nnz( full(E) - E_expected )
%    fprintf(2, 'Unable to verify results\n');
%else
%    disp(E);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph Challenge benchmark
% Developer : Dr. Siddharth Samsi (ssamsi@mit.edu)
%
% MIT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) <2017> Massachusetts Institute of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

