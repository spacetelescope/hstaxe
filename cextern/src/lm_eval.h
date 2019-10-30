/**
 */
typedef struct {
  double* user_t;
  double* user_y;
  double (*user_func)( const double user_t_point, const double* par );
} lm_data_type;

typedef struct {
  double* user_t;
  double* user_y;
  double* fpars;
  double (*user_func)(const double user_t_point, const double* par,const  double* fpar);
} lm_data_fpar_type;

void lm_evaluate_default( double* par, int m_dat, double* fvec,
                          void *data, int *info );

void lm_evaluate_fpar( double* par, int m_dat, double* fvec,
		       void *data, int *info );

void lm_print_default( int n_par, double* par, int m_dat, double* fvec,
                       void *data, int iflag, int iter, int nfev );
void lm_print_nothing( int n_par, double* par, int m_dat, double* fvec,
                       void *data, int iflag, int iter, int nfev );
