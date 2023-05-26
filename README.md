# M-exdex
[![View M-exdex (Estimation of the Extremal Index) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/130094-m-exdex-estimation-of-the-extremal-index)
 
 Matlab/Octave distribution of the exdex R package (Estimation of the Extremal Index):
 
 
 https://github.com/paulnorthrop/exdex/tree/master
 
 https://paulnorthrop.github.io/exdex/
 
 
 All code has only been minimally tested (but results mirror that of the R package) - use with caution
 
 
At present, only the functions for the extremal index calculation have been completed, in the future I will look to translating the rest of the package, such as the diagnostic plots. For full details about the theory, functions, and inputs, please refer to the original package.
 
 
Syntax is exdex.function(), the currently working functions are:

```matlab:Code
theta = exdex.spm(data, b, bias_adjust, constrain, varN, which_dj, nv);
theta = exdex.kgaps(data, u, k, inc_cens, nv);
theta = exdex.dgaps(data, u, D, inc_cens, nv);
theta = exdex.iwls(data, u, maxit, nv);
% 'nv' is name value argument 'disp' to specify whether to display output to command window or not.
% default is to display, specify 'disp', 'n' to not display.
```

All credit and thanks to authors of the original R package.
