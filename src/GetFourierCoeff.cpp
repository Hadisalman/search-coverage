#include "GetFourierCoeff.h"

void GetFourierCoeff(opt1 *opt,MatrixXd X,MatrixXd Y)
{
	int a;
	float PI=3.14159265359;
	int Nkx=opt->erg.Nkx;
	int Nky=opt->erg.Nky;
	MatrixXd muk(Nkx,Nky);
	MatrixXd temp1(1 , Nkx);
	MatrixXd temp2(1 , Nky);

	temp1<<1/sqrt(0.5),MatrixXd::Ones(1,Nkx-1);
	temp2<<1/sqrt(0.5),MatrixXd::Ones(1,Nky-1);
	
	//Initializing Fourier coefficents of prior
	
	muk=MatrixXd::Zero(Nkx,Nky);

	opt->erg.HK=sqrt(opt->Lx*(opt->Ly))*((sqrt(0.5)*(temp2).transpose())*(sqrt(0.5)*(temp1) ) );
	
	cout<<(opt->erg.HK).size()<<endl;
 	
 	for(int kx=0;kx< Nkx;kx++)
 	{	
 		for(int ky=0;ky<Nky;ky++)
 		{	a=opt->erg.HK(ky,kx);
 			muk(kx,ky) = ((((opt->erg.mu).cwiseProduct((PI*(kx)*(X)/opt->Lx).array().cos().matrix()))).cwiseProduct(((PI*(ky)*(Y)/opt->Ly).array().cos().matrix())).sum())/a;
 			
 		}
 	}
 	
 	opt->erg.muk=muk;

}