#include "initialization.h"

void initialization(opt1 *opt,pos *pose)
{	
	float PI;
	PI=3.141592;
	int i;
	//Domain bounds
	opt->dbxmin=0;
	opt->dbxmax=300;
	opt->dbymin=0;
	opt->dbymax=300;
	
	opt->Lx=opt->dbxmax-opt->dbxmin;
	opt->Ly=opt->dbymax-opt->dbymin;
	
	//Number of agents				
	opt->nagents=3;

	//Velocity boundes(lower and upper)
	opt->vlb=0.1*MatrixXd::Ones(opt->nagents,1);
	opt->vub=5*MatrixXd::Ones(opt->nagents,1);

	opt->wlb=-0.5*MatrixXd::Ones(opt->nagents,1);
	opt->wub=0.5*MatrixXd::Ones(opt->nagents,1);

	//initializing agent locations
	pose->x<<100,200,50;
	pose->y<<50,75,250;
	pose->theta<<(PI/4),(PI/4),(PI/4);

	//ergodicity parameters

	opt->erg.s=1.5;  //%parameter of soblev norm

	opt->erg.Nkx=50;
	opt->erg.Nky=50;

	VectorXd v(opt->erg.Nkx);
	MatrixXd m;

	m=MatrixXd::Ones(1,opt->erg.Nkx);

	for(i=0;i<opt->erg.Nkx;i++)     
	{
		v(i) = i;

	}

	opt->erg.KX= v*m;
	opt->erg.KY=MatrixXd::Ones(opt->erg.Nkx,1)*(v.transpose());
	opt->erg.LK=((MatrixXd::Ones(opt->erg.Nky,opt->erg.Nky)+ (opt->erg.KX).array().pow(2).matrix() + (opt->erg.KY).array().pow(2).matrix()).array().pow(opt->erg.s).matrix()).cwiseInverse();

	// simulation params
	opt->simNsteps=5000;
	opt->simdt=(1e-1);

}