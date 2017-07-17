#include "SMC_Update.h"

void SMC_Update(pos *pose, MatrixXd *Ck, float time,opt1 opt)
{
	float PI;
	PI=3.141592;
	float xrel;
	float yrel;
	float Bjx;
	float Bjy;
	float GammaV;
	float GammaW;
	float v;
	float w;
	MatrixXd HK;
	MatrixXd temp1;
	MatrixXd temp2;
	MatrixXd temp3;
	


	HK=opt.erg.HK;
	
	for(int i=0;i< opt.nagents;i++)
	{
		xrel=pose->x(i)-opt.dbxmin;
		yrel=pose->y(i)-opt.dbymin;

		temp1=((opt.erg.KX)*(PI*xrel/opt.Lx)).array().cos().matrix();
		temp2=((opt.erg.KY)*(PI*yrel/opt.Ly)).array().cos().matrix();
		
		(*Ck)=(*Ck) +(((temp1).cwiseProduct((temp2)*(opt.simdt))).cwiseQuotient((opt.erg.HK).transpose()));
		
	}

	for (int i = 0; i < opt.nagents; ++i)
	{
		xrel=pose->x(i)-opt.dbxmin;
		yrel=pose->y(i)-opt.dbymin;
		temp1=(opt.erg.LK).cwiseQuotient((opt.erg.HK).transpose());
		temp2=(*Ck)-opt.nagents*(time)*(opt.erg.muk);
		temp3=(opt.erg.KX*(PI*xrel/opt.Lx)).array().sin().matrix() ;
		
		Bjx=( ( (( temp1 ).cwiseProduct(temp2) ).cwiseProduct((-opt.erg.KX*(PI/opt.Lx)).cwiseProduct(temp3) )).cwiseProduct(((opt.erg.KY)*(PI*yrel/opt.Ly)).array().cos().matrix())).sum();
		
		Bjy=( ( (( (opt.erg.LK).cwiseQuotient((opt.erg.HK).transpose() ) ).cwiseProduct((*Ck)-opt.nagents*time*(opt.erg.muk)) ).cwiseProduct((-opt.erg.KY*(PI/opt.Ly)).cwiseProduct((opt.erg.KX*(PI*xrel/opt.Lx)).array().cos().matrix()) )).cwiseProduct((opt.erg.KY*(PI*yrel/opt.Ly)).array().sin().matrix())).sum();
		
		GammaV=Bjx*cos(pose->theta(i))+Bjy*sin(pose->theta(i));
		GammaW=-Bjx*sin(pose->theta(i))+Bjy*cos(pose->theta(i));
		
		//Updating agent location based on SMC feedback control law
		if (GammaV>=0)
			v=opt.vlb(i);
		else 
			v=opt.vub(i);

		if (GammaW>=0)
			w=opt.wlb(i);
		else 
			w=opt.wub(i);

		//velocity motion model
		if (abs(w)<1e-10)
		{
			pose->x(i)=	pose->x(i) +v*(opt.simdt)*cos(pose->theta(i));
			pose->y(i)=	pose->y(i) +v*(opt.simdt)*sin(pose->theta(i));
		}
		else
		{
			pose->x(i)=	pose->x(i) +(v/w)*(sin(pose->theta(i)+w*(opt.simdt)) - sin(pose->theta(i)));
			pose->y(i)=	pose->y(i) +(v/w)*(cos(pose->theta(i)) - cos(pose->theta(i)+w*(opt.simdt)));
		}

		pose->theta(i)=pose->theta(i) +w*(opt.simdt);

	}
}