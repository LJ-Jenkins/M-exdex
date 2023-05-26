# M-exdex
 Matlab/Octave distribution of the exdex R package:
 
 
 https://github.com/paulnorthrop/exdex/tree/master
 
 https://paulnorthrop.github.io/exdex/
 
 
 All code has only been minimally tested - use with caution -
 
 
At present, only the functions extremal index calculation functions have been completed, in the future I will look to translating the rest of the package. For full details about the theory, functions, and inputs, please refer to the original package.
 
 
Syntax is exdex.function(), currently working functions:

```matlab:Code
theta = exdex.spm(data, b, bias_adjust, constrain, varN, which_dj, nv);
theta = exdex.kgaps(data, u, k, inc_cens, nv);
theta = exdex.dgaps(data, u, D, inc_cens, nv);
theta = exdex.iwls(data, u, maxit, nv);
% 'nv' is name value argument 'disp' to specify whether to display output to command window or not.
% default is to display, specify 'disp', 'n' to not display.
```
