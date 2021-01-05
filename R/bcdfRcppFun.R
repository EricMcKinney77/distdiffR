cppFunction('
  NumericVector bcdf(NumericMatrix data, NumericMatrix eval) {
    NumericVector x(eval.nrow());
    int i, j;
    for(i = 0; i < eval.nrow(); ++i) {
      x(i) = 0;
      for(j = 0; j < data.nrow(); ++j) {
        if( data(j, 0) <= eval(i, 0) && data(j, 1) <= eval(i, 1) ) x(i) = x(i) + 1.0;
      }
      x(i) = x(i) / data.nrow();
    }
    return x;
  }
')
