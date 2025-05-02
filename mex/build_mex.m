target = 'D:\Trenager\SEUIF97\mex\mex_seuif97.c';
includepath = '-ID:\Trenager\SEUIF97\src\common';
libpath = 'D:\Trenager\SEUIF97\build\Release\seuif97.lib';
%mex('-v','-g','-R2018a',target,includepath,libpath);
mex('-v','-R2018a',target,includepath,libpath);
clear