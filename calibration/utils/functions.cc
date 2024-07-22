#include <utils/functions.h>

void squeeze(array& arr) {
    arr.resize(arr.size(), 1);
};

array rowToMatrix(const array& arr, unsigned n) {
    return arr.replicate(n, 1);
};

array columnToMatrix(const array& arr, unsigned m) {
    return arr.replicate(1, m);
};