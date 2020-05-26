os = computer();
install_path = [userpath,'/SNPLIB'];
switch os
    case 'PCWIN64'
        mkl_root = 'C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl';
    case 'GLNXA64'
        mkl_root = '/opt/intel/mkl';
        intel_root = '/opt/intel';
    case 'MACI64'
        mkl_root = '/opt/intel/mkl';
        intel_root = '/opt/intel';
end
mkl_include = ['-I"',mkl_root,'/include"'];
switch os
    case 'PCWIN64'
        mkl_seq_libs = ['LINKLIBS=$LINKLIBS ','/LIBPATH:"',mkl_root,'\lib\intel64" ' ,'mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib '];
        cxxoptim = '';
    case 'GLNXA64'
        mkl_seq_libs = ['LINKLIBS=$LINKLIBS -Wl,--start-group ',mkl_root,'/lib/intel64/libmkl_intel_lp64.a ',mkl_root,'/lib/intel64/libmkl_sequential.a ',mkl_root,'/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl'];
        cxxoptim = 'CXXOPTIMFLAGS=$CXXOPTIMFLAGS -O3 -std=c++11 -march=native -pipe -fPIC -flto';
    case 'MACI64'
        mkl_seq_libs = ['LINKLIBS=$LINKLIBS ',mkl_root,'/lib/libmkl_intel_lp64.a ',mkl_root,'/lib/libmkl_sequential.a ',mkl_root,'/lib/libmkl_core.a',' -lpthread -lm -ldl'];
        cxxoptim = 'CXXOPTIMFLAGS=$CXXOPTIMFLAGS -O3 -std=c++11 -march=native -pipe -fPIC -flto -mmacosx-version-min=10.7';
end
mkdir(install_path);
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcAdjustedAF_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcAdjustedMAF_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcAdjustedGRM_.cc','src/adjusted_grm.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcAdmixedGRM_.cc','src/adjusted_grm.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcAlleleFrequencies_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcCCAGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcCCAReplication_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'matlab/mexfiles/CalcGCTADiagonal_.cc','src/grm.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcGRMMatrix_.cc','src/grm.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcUGRMMatrix_.cc','src/ugrm.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcIBSConnection_.cc','src/ibs.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcIBSMatrix_.cc','src/ibs.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcKINGMatrix_.cc','src/king.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/FindUnrelatedGroup_.cc','src/king.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcLinearRegressionGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcLogisticGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcMissing_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcMultiLMM_REML_.cc','src/genetic_variances.cc','src/snp.cc','src/uni_lmm.cc','src/multi_lmm.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcMultiLMM_RML_.cc','src/genetic_variances.cc','src/snp.cc','src/uni_lmm.cc','src/multi_lmm.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcUniLMMGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_MKL',mkl_include,mkl_seq_libs,'matlab/mexfiles/CalcUniLMM_.cc','src/genetic_variances.cc','src/snp.cc','src/uni_lmm.cc','src/multi_lmm.cc');
mex(cxxoptim,'matlab/mexfiles/GenerateAdmixedIndividuals_.cc','src/simulations.cc');
mex(cxxoptim,'matlab/mexfiles/GenerateIndividuals_.cc','src/simulations.cc');
mex(cxxoptim,'matlab/mexfiles/GeneratePairwiseSiblings_.cc','src/simulations.cc');
mex(cxxoptim,'matlab/mexfiles/UpdateAf_.cc','src/simulations.cc');
mex(cxxoptim,'matlab/mexfiles/UnpackGeno_.cc','src/data_manage.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/UnpackGRMGeno_.cc','src/data_manage.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/UnpackUGeno_.cc','src/data_manage.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/FlipGeno_.cc','src/data_manage.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/Keep_.cc','src/data_manage.cc','src/snp.cc');
copyfile('matlab/@SNPLIB/*.m',[install_path,'/@SNPLIB']);
copyfile('matlab/GLMM/*.m',[install_path,'/GLMM']);
movefile(['*.',mexext],[install_path,'/mexfiles']);
cd(install_path);
clear all;