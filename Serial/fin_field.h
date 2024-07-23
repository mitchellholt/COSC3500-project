// Compute a*b mod p
int mulmod(int a, int b, int p);

// Compute a + b mod p
int addmod(int a, int b, int p);

// Compute a - b mod p
int submod(int a, int b, int p);


// Return a primitive nth root of unity in Z/p if one exists. If not, return 0
int primitive_root(int n, int p);

// Find the (multiplicative) inverse of n in Z/p
int invmod(int n, int p);
