#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;
using namespace std;

typedef struct ergg
{
	float s;
	int Nkx;
	int Nky;
	MatrixXd LK;
	MatrixXd KX;
	MatrixXd KY;
 	MatrixXd mu;
 	MatrixXd muk;
 	MatrixXd HK;

}ergg;


typedef struct opt1

{ 
	int dbxmin;
	int dbxmax;
	int dbymin;
	int dbymax;
	int Lx;
	int Ly;
	int nagents;
	ergg erg;
	int simNsteps;
	float simdt;
	Vector3d vlb;
	Vector3d vub;
	Vector3d wlb;
	Vector3d wub;
	
}opt1;

typedef struct pos
{
	Vector3d x;
	Vector3d y;
	Vector3d theta;
	
}pos;

#endif
