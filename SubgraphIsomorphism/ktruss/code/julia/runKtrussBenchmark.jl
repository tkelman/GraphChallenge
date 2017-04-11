
include("ktruss.jl");

inc_mtx_file = "../../../data/amazon0302_inc.tsv"

E_expected =  [1  1  0  0  0
               0  1  1  0  0
               1  0  0  1  0
               0  0  1  1  0
               1  0  1  0  0
               0  0  0  0  0];


@time E = ktruss(inc_mtx_file, 3);

#if sum( E - E_expected ) > 0
#    println("Unable to verify results");
#else
#    println("passed");
#    println(E);
#end

writedlm("amazon0302_result.tsv", hcat(findnz(E)...), '\t')

#######################################################
# Graph Challenge Benchmark
# Developer : Dr. Siddharth Samsi (ssamsi@mit.edu)
#
# MIT
#########################################################
# (c) <2017> Massachusetts Institute of Technology
########################################################


