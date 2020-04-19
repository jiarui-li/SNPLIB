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
system(['cmake.exe ./matlab -DUSE_MKL=ON -DMKL_ROOT="',mkl_root,'"']);
system('cmake.exe --build . --config Release');

copyfile('matlab/@SNPLIB/*.m',[install_path,'/@SNPLIB']);
copyfile('matlab/GMLM/*.m',[install_path,'/GMLM']);
copyfile('matlab/misc/*.m',[install_path,'/misc']);
movefile(['*.',mexext],[install_path,'/mexfiles']);
cd(install_path);
clear all;