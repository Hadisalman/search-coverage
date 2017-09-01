#include <iostream>
#include <ctime> 
#include <Eigen/Dense>  
#include "structures.h"
#include "initialization.h"
#include "GetFourierCoeff.h"
#include "GenerateUtilityMap.h"
#include "calculate_Ergodicity.h"
#include "SMC_Update.h"
#include "ros/ros.h"
#include "geometry_msgs/Pose2D.h"  //we can change Int8 by any std::msgs types
#include "std_msgs/Float32.h" 
#include "std_msgs/Float32MultiArray.h" 

#include <sstream>

using namespace Eigen;
using namespace std;

const float PI = 3.1415927;
//#define PI 3.14159265;



int main(int argc, char **argv) 
{	
	
//////////////////////// ROS inititazaiton*/
	ros::init(argc, argv, "ergodic_search_node");

	ros::NodeHandle n;
 
    ros::Publisher pose_pub = n.advertise<std_msgs::Float32MultiArray>("robot_traj", 1000);
    ros::Publisher iterations_pub = n.advertise<std_msgs::Float32>("iterations_num", 1000);
  	ros::Rate loop_rate(1000);
//////////////////////////////////////////////

	time_t tstart;
	time_t tend;
	pos pose;
	opt1 opt;
	int addnoise=0;

	initialization(&opt,&pose);

	MatrixXd X(opt.dbymax,opt.dbxmax);
	MatrixXd Y(opt.dbymax,opt.dbxmax);
	MatrixXd Ck;
	MatrixXd Z;
	MatrixXd Ergodicity_Metric(opt.simNsteps,1);
	MatrixXd xposition(opt.simNsteps,opt.nagents);
	MatrixXd yposition(opt.simNsteps,opt.nagents);
	MatrixXd ck;
	float t;
	int Nsteps=opt.simNsteps;
 	float dt=opt.simdt;
 	int size;
 	size=X.size();
 	VectorXd G(size);

	
	GenerateUtilityMap(&X,&Y,&G,opt,addnoise);
 	
 	Z=MatrixXd::Zero(opt.dbymax-1,opt.dbxmax-1);

 	Map<MatrixXd> M2(G.data(), opt.dbymax,opt.dbxmax);

 	opt.erg.mu=M2;
 	opt.erg.mu=opt.erg.mu/((opt.erg.mu).sum());
 	
	GetFourierCoeff(&opt,X,Y);

 	//run simulation
 	
 	Ck=MatrixXd::Zero(opt.erg.Nkx,opt.erg.Nky);

 	Ergodicity_Metric=MatrixXd::Zero(Nsteps,1);
 	xposition=MatrixXd::Zero(Nsteps,opt.nagents);
 	xposition=MatrixXd::Zero(Nsteps,opt.nagents);
 	
	std_msgs::Float32MultiArray pose_msg;
 	std_msgs::Float32 iterations_msg;
 	

 	tstart = time(0);
	for (int i = 0; i < Nsteps; ++i)
	 	{	

		iterations_msg.data=i;
	 	iterations_pub.publish(iterations_msg);
 		
 		t=(i+1)*dt;
 		SMC_Update(&pose,&Ck,t,opt);
 		
 		for(int j=0;j<opt.nagents;j++)
 		{	
 			xposition(i,j)=pose.x(j);
 			yposition(i,j)=pose.y(j);
 		}
 		// pose_msg.x =  pose.x(2);
 		// pose_msg.y =  pose.y(2);
 		if(i%20==0){
 		pose_msg.data.clear();
 		// cout<<pose.x(0)<<endl;
 		pose_msg.data.push_back(pose.x(0));
 		pose_msg.data.push_back(pose.y(0));
 		pose_msg.data.push_back(pose.x(1));
 		pose_msg.data.push_back(pose.y(1));
 		pose_msg.data.push_back(pose.x(2));
 		pose_msg.data.push_back(pose.y(2));
 		// pose_pub.publish(pose_msg);
 		pose_pub.publish(pose_msg);
 		}

 		ck=Ck/opt.nagents/t;

 		Ergodicity_Metric(i,0) = calculate_Ergodicity(ck,opt);
	 	
	 		
	 	}
	 	
	 tend = time(0); 
	cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl; 	

  return 0;
}