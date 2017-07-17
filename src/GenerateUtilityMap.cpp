
#include "GenerateUtilityMap.h"

void GenerateUtilityMap(MatrixXd *X,MatrixXd *Y,VectorXd *G,opt1 opt,int addnoise)
{
	int j=0;
	int xdel;
	int ydel;
	VectorXd xRange(300);
	VectorXd yRange(300);
	float PI;
	PI=3.141592;
	VectorXd m1(2);
	VectorXd m2(2);
	VectorXd m3(2);
	MatrixXd s1;
	MatrixXd s2;
	MatrixXd s3;
	int size;
	size=(*X).size();
	MatrixXd x(size,2);	
	
	VectorXd G1(size);
	VectorXd G2(size);
	VectorXd G3(size);

	xdel=1;
	ydel=1;
		
	for(int i = opt.dbxmin; i<=opt.dbxmax-xdel; i+=xdel)
	{
		xRange(j)=i;
		j++;
	}
	
	j=0;
	
	for(int i = opt.dbymin; i<=opt.dbymax-xdel; i+=xdel)
	{
		yRange(j)=i;
		j++;
	}
	
	//meshgrid

	for(int i = 0; i<=opt.dbymax-ydel; i++)
	{
		(*X).row(i)=xRange.transpose();

	}

	for(int i = 0; i<=opt.dbxmax-xdel; i++)
	{
		(*Y).col(i)=yRange;
	}

	m1<<100,200;
	s1=150*MatrixXd::Identity(2,2);

	m2<<220, 200;
	s2=800*MatrixXd::Identity(2,2);

	m3<<120,120;
	s3=600*MatrixXd::Identity(2,2);

	VectorXd temp1(Map<VectorXd>((*X).data(), size));
	VectorXd temp2(Map<VectorXd>((*Y).data(), size));
	x<<temp1,temp2;

	for (int i = 0; i < size; ++i)
	{

		G1(i)=(1/(PI*sqrt(s1.determinant())))*(exp(-0.5*(x.row(i)-m1.transpose())*(s1.inverse())*((x.row(i)-m1.transpose()).transpose())));	
		G2(i)=(1/(PI*sqrt(s2.determinant())))*(exp(-0.5*(x.row(i)-m2.transpose())*(s2.inverse())*((x.row(i)-m2.transpose()).transpose())));	
		G3(i)=(1/(PI*sqrt(s3.determinant())))*(exp(-0.5*(x.row(i)-m3.transpose())*(s3.inverse())*((x.row(i)-m3.transpose()).transpose())));	
		
	    
	}	

     (*G)=G1+3*G2+G3;
     (*G)=(*G).cwiseMax(0);
     (*G)=(*G)/((*G).maxCoeff());
}