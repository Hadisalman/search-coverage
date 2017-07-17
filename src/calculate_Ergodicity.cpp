
#include "calculate_Ergodicity.h"

using namespace Eigen;

double calculate_Ergodicity(MatrixXd ck, opt1 opt)
{
	double Ergodicity_Metric;
	
	Ergodicity_Metric=((opt.erg.LK).cwiseProduct(((ck-opt.erg.muk)).array().pow(2).matrix())).sum();

	return Ergodicity_Metric;
}
