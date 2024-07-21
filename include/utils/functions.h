#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <utils/__types__.h>

// not a traditional (i.e., NumPy's) squeeze since it doesnt collapse
// a dimension, but it creates an [equivalent N x 1 column Eigen array] ~
// not a traditional (i.e., NumPy's) squeeze since it doesnt collapse
// a dimension, but it creates an [equivalent N x 1 column Eigen array] ~
void squeeze(array &arr);

array rowToMatrix(const array &arr, unsigned n = 15);

array colToMatrix(const array &arr, unsigned m = 11);

#endif  // FUNCTIONS_H