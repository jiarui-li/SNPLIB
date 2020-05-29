os = computer();
install_path = [userpath,'/SNPLIB'];
switch os
    case 'PCWIN64'
        openblas_root = 'C:\opt';
        openblas_include = ['-I"',openblas_root,'\include\openblas"'];
        clang_library_path = 'C:\ProgramData\Miniconda3\Library\lib';
        openblas_lib = ['LINKLIBS=$LINKLIBS ','/LIBPATH:"',openblas_root,'\lib" ','openblas.lib /LIBPATH:"',clang_library_path,'" flangmain.lib flang.lib flangrti.lib ompstub.lib'];
        cxxoptim= '';
    case 'GLNX64'
        openblas_root = '/opt/OpenBLAS';
        openblas_include = ['-I"',openblas_root,'\include\"'];
        openblas_lib = ['LINKLIBS=$LINKLIBS -L"',openblas_root,'/lib" -lopenblas -lgfortran'];
        cxxoptim = 'CXXOPTIMFLAGS=$CXXOPTIMFLAGS -O3 -std=c++11 -march=native -pipe -fPIC -flto -DUSE_OPENBLAS';
    case 'MACI64'
        openblas_root = '/opt/OpenBLAS';
        openblas_include = ['-I"',openblas_root,'\include\"'];
        openblas_lib = ['LINKLIBS=$LINKLIBS -L"',openblas_root,'/lib" -lopenblas -lgfortran'];
        cxxoptim = 'CXXOPTIMFLAGS=$CXXOPTIMFLAGS -O3 -std=c++11 -march=native -pipe -fPIC -flto -DUSE_OPENBLAS';
end
mkdir(install_path);
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcAdjustedAF_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcAdjustedMAF_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcAdjustedGRM_.cc','src/adjusted_grm.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcAdmixedGRM_.cc','src/adjusted_grm.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcAlleleFrequencies_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcCCAGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcCCAReplication_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'matlab/mexfiles/CalcGCTADiagonal_.cc','src/grm.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcGRMMatrix_.cc','src/grm.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcUGRMMatrix_.cc','src/ugrm.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcIBSConnection_.cc','src/ibs.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcIBSMatrix_.cc','src/ibs.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/CalcKINGMatrix_.cc','src/king.cc','src/snp.cc');
mex(cxxoptim,'matlab/mexfiles/FindUnrelatedGroup_.cc','src/king.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcLinearRegressionGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcLogisticGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcMissing_.cc','src/statistics.cc','src/snp.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcMultiLMM_REML_.cc','src/genetic_variances.cc','src/snp.cc','src/uni_lmm.cc','src/multi_lmm.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcMultiLMM_RML_.cc','src/genetic_variances.cc','src/snp.cc','src/uni_lmm.cc','src/multi_lmm.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcUniLMMGWAS_.cc','src/gwas.cc','src/snp.cc','src/uni_lmm.cc');
mex(cxxoptim,'-DUSE_OPENBLAS',openblas_include,openblas_lib,'matlab/mexfiles/CalcUniLMM_.cc','src/genetic_variances.cc','src/snp.cc','src/uni_lmm.cc','src/multi_lmm.cc');
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
movefile(['*.',mexext],[install_path,'/mexfiles']);
cd(install_path);
clear all;