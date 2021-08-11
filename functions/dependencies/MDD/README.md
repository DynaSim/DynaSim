# MDD
MATLAB Multidimensional Dictionary class

**Summary**: MDD allows the implementation of multidimensional Python-esque dictionaries in MATLAB.

**Description**: MDD can be interpreted in several different ways. 
- A map/dictionary that associates multiple keys with a single value
- An N-dimensional table (a table is equivalent to an MDD object in 2-dimensions)
- A matrix or cell array that can be indexed by using strings and regular expressions

**WARNING**: This is a work in progress. The core code is implemented, but documentation is missing. I will be updating it more over the coming months.

This code is similar to the multidimensional map function implemented by David Young. However, this implementation does not use MATLAB maps and instead is adds functionality to traditional MATLAB matrices and cell arrays.

**Related commands**: 
- MATLAB Map Containers (https://www.mathworks.com/help/matlab/map-containers.html)
- Multidimensional implementation of this by David Young (https://www.mathworks.com/matlabcentral/fileexchange/33068-a-multidimensional-map-class)
