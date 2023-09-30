extern "C"
{
  //http://www.netlib.org/lapack/double/dgesv.f
  void dgesv_(int* n, int *nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
}
 
int main()
{
  int n = 2;
  int nrhs = 2;
  double a[] = {1,3,2,4}; //colum-major mode 
  int lda = n;
  int* ipiv = new int[n];
  double b[] = {19, 43, 22, 50}; //column-major mode
  int ldb = n;
  int info;
   
  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
 
  delete [] ipiv;
  return 0;
}