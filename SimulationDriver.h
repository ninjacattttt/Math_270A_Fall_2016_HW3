#ifndef _SIMULATION_DRIVER_H_
#define _SIMULATION_DRIVER_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <fstream>
#include "LagrangianForce.h"

namespace FILE_IO{
  inline void Write_Binary(std::string directory, std::string name,Eigen::VectorXd& v){
    std::string filename(directory+std::string("/")+name+std::string(".binary"));
    std::ofstream outdata(filename.c_str(),std::ios::out|std::ios::binary);
    for(int i=0;i<v.size();i++)
      outdata.write(reinterpret_cast<char*>(&v(i)),sizeof(double));
    outdata.close();
  }

  inline void Read_Binary(std::string directory, std::string name,Eigen::VectorXd& v){
    std::string filename(directory+std::string("/")+name+std::string(".binary"));
    std::ifstream indata(filename.c_str(),std::ios::in|std::ios::binary);
    for(int i=0;i<v.size();i++)
      indata.read(reinterpret_cast<char*>(&v(i)),sizeof(double));
    indata.close();
  }

  inline void Write_Binary(std::string directory, std::string name,Eigen::VectorXf& v){
    std::string filename(directory+std::string("/")+name+std::string(".binary"));
    std::ofstream outdata(filename.c_str(),std::ios::out|std::ios::binary);
    for(int i=0;i<v.size();i++)
      outdata.write(reinterpret_cast<char*>(&v(i)),sizeof(float));
    outdata.close();
  }

  inline void Read_Binary(std::string directory, std::string name,Eigen::VectorXf& v){
    std::string filename(directory+std::string("/")+name+std::string(".binary"));
    std::ifstream indata(filename.c_str(),std::ios::in|std::ios::binary);
    for(int i=0;i<v.size();i++)
      indata.read(reinterpret_cast<char*>(&v(i)),sizeof(float));
    indata.close();
  }

  inline void Write_DAT_File(std::string file,const Eigen::VectorXf& array){
    FILE* fpointer;
    fpointer=fopen(file.c_str(),"w");
    for(int i=0;i<array.size();i++)
        fprintf(fpointer,"%f\n",array(i));
    fclose(fpointer);
}

inline void Write_DAT_File(std::string file,const Eigen::VectorXd& array){
    FILE* fpointer;
    fpointer=fopen(file.c_str(),"w");
    for(int i=0;i<array.size();i++)
        fprintf(fpointer,"%g\n",array(i));
    fclose(fpointer);
}
}

namespace JIXIE{
template <class T>
class SimulationParameters{
public:
	T final_time;
	int frames_per_second;
	T dt;
	std::string output_dir;
	SimulationParameters(){}		
};
	
template <class T>
class SimulationDriver{
public:
	T time,dt,dt_target,dt_frame,final_time,dt_min;					
	int frames_per_second,current_frame;
	std::string output_directory;

	SimulationDriver(const T final_time_input,const int frames_per_second_input,const T dt_input,std::string& output_dir):final_time(final_time_input),frames_per_second(frames_per_second_input),dt(dt_input),dt_target(dt),dt_min((T)1e-6),output_directory(output_dir){
    dt_frame=(T)1/(T)frames_per_second;							
    if(dt_target > dt_frame)									//dt_target is the largest possible dt to ensure stability
      	dt_target=dt_frame;										//if dt_frame is smaller then use as dt_target and dt
    dt=dt_target;
  	}

  	SimulationDriver(SimulationParameters<T>& parameters):SimulationDriver(parameters.final_time,parameters.frames_per_second,parameters.dt,parameters.output_dir){}

  	virtual void Set_Dt(bool& write_frame){
    	dt=dt_target;
    	if((T)(current_frame + 1) * dt_frame - time < dt_min + dt){  		//If next next step will be smaller than dt_min then take one bigger timestep
	  		dt=(T)(current_frame + 1) * dt_frame - time;
	  		current_frame++;												//jump to the next frame
      		write_frame=true;												//write frame
    	}		
    	else if(time + dt > (T)(current_frame + 1) * dt_frame){				//If next step jumps over frame time then take one smaller timestep 
      		dt=(T)(current_frame + 1) * dt_frame - time;
      		current_frame++;
      		write_frame=true;
    	}
    	else write_frame=false;												//If not then no need to write frame
  	}

  	virtual void Initialize(){
    	time=(T)0;
    	WriteState(current_frame);
  	}

  	void WriteState(const int number){
    	std::ofstream outdata;
    	std::string simulation_data_filename(output_directory+std::string("/simulation_info.dat"));
    	outdata.open(simulation_data_filename.c_str());
    	outdata << number << std::endl;
    	outdata.close();
    	Write_State(number,simulation_data_filename);
  	}
	
	virtual void WriteObj(const int number){}
	virtual void InitObj(){}

  	virtual void Write_State(const int number,std::string& simulation_data_filename){}
  	virtual bool Read_State(const int number,std::string& simulation_data_filename){return false;}

  	virtual void Advance_One_Time_Step(const bool verbose){
    	time+=dt;
  	}

  	void RunSimulation(const bool verbose=false){
		Initialize();
		InitObj();
    	while(time<final_time){
      	  	if(verbose)
        		std::cout << "Time = " << time << ", frame = " << current_frame << ", dt = " << dt << std::endl;
      	  	bool write_to_file=false;
      	  	Set_Dt(write_to_file);
      	  	Advance_One_Time_Step(verbose);
      	  	if(write_to_file){
				WriteState(current_frame);
				WriteObj(current_frame);
			}
				
    	}
  	}
};

template <class T>
class ElasticityParameters:public SimulationParameters<T>{
public:
	int N;
	T a;
	T dX;
	T rho;
	T k;
	T Newton_tol;
	int max_newton_it;
	
	ElasticityParameters(){}
};

template <class T>
class ElasticityDriver: public SimulationDriver<T>{
	using SimulationDriver<T>::output_directory;
	using SimulationDriver<T>::time;
  	using SimulationDriver<T>::dt;
  	typedef Eigen::Matrix<T,Eigen::Dynamic, 1> TVect;
  	int N;
  	T a,dX;
  	T rho,k;
  	TVect x_n,x_np1,v_n,x_hat,residual,mass,delta;
  	T Newton_tol;
  	int max_newton_it;
  	ConstitutiveModel<T>* cons_model;
  	LagrangianForces<T>* lf;
  	SymmetricTridiagonal<T> be_matrix;
public:
  /*ElasticityDriver(const T final_time_input,const int frames_per_second_input,const T dt_input,const int N_input,const T a_input,const T dX_input,std::string& output_dir):
  SimulationDriver<T>(final_time_input,frames_per_second_input,dt_input,output_dir),N(N_input),a(a_input),dX(dX_input),
  rho((T)1),k((T)100),x_n(N_input),x_np1(N_input),v_n(N_input),x_hat(N_input),residual(N_input),mass(N_input),delta(N_input),
  Newton_tol((T)1e-5),max_newton_it(10),be_matrix(N_input){
    cons_model=new LinearElasticity<T>(k);*/
	
	ElasticityDriver(ElasticityParameters<T>& parameters):SimulationDriver<T>(parameters),N(parameters.N),a(parameters.a),dX(parameters.dX),rho(parameters.rho),k(parameters.k),x_n(parameters.N),x_np1(parameters.N),v_n(parameters.N),x_hat(parameters.N),residual(parameters.N),mass(parameters.N),delta(parameters.N),Newton_tol(parameters.Newton_tol),max_newton_it(parameters.max_newton_it),be_matrix(parameters.N){
		cons_model=new NeoHookean<T>(k);								
		//cons_model=new LinearElasticity<T>(k);
    	lf=new FEMHyperelasticity<T>(a,dX,N,*cons_model);
	}

  	~ElasticityDriver(){
    	delete cons_model;
    	delete lf;
  	}

  	void Initialize(){													//set initial positions and velocity	
  		for(int i=0;i<N;i++){
			T x=(a+(T)i*dX);
			x_n(i)=(T).7*x;
      		v_n(i)=(T)0;
    	}
    	//intialize mass lumped mass matrix from density
		
		mass(0)=(T)-0.5*rho*dX;											//there's no N_0 basis function 
		mass(1)=(T)-(1/6)*rho*dX;											
	
    	for(int e=0;e<N-1;e++){													
      		mass(e)+=(T).5*rho*dX;
      		mass(e+1)+=(T).5*rho*dX;
  		}
    	SimulationDriver<T>::Initialize();
  	}

  	virtual void Advance_One_Time_Step(const bool verbose){					
    	time+=dt;
    	x_hat=x_n+dt*v_n;
    	x_np1=x_n;//initial guess

    	for(int it=1;it<max_newton_it;it++){
      	  	residual=mass.asDiagonal()*(x_hat-x_np1);
    	  	lf->AddForce(residual,x_np1,dt*dt);							//residual=-g(x_np1)
    	  	
			residual(0)=0;												//mass(0)=0, f(0) should be 0
		  
		  	T norm=(T)0;
		  	for(int i=1;i<N;i++) 
			  	norm+=residual(i)*residual(i)/mass(i);
    	  	norm=sqrt(norm);
	  
		  	if(verbose)
      		  	std::cout << "Newton residual at iteration " << it << " = " << norm << std::endl;
		  	if(norm<Newton_tol){
			  	Exit_BE();
			  	return;
		 	}
	  
    	  	be_matrix.SetToZero();										//be_matrix=M-dt*dt*df(x_np1)
    	  	for(int i=0;i<N;i++) 
				be_matrix(i,i)=mass(i);
    	  	lf->AddForceDerivative(be_matrix,x_np1,-dt*dt);
		  
		  	be_matrix(0,0)=1;											//ensure delta(0)=0 so that x(0)=0
		  	be_matrix(0,1)=0;
		  	be_matrix(1,0)=0;
		  
    	  	//be_matrix.QRSolve(delta,residual);
			
			be_matrix.GMRES(delta,residual);							//GMRES method
			
    	  	x_np1+=delta;
		}
		
		Exit_BE();
	}
  
  	void Exit_BE(){
  		v_n=((T)1/dt)*(x_np1-x_n);
		x_n=x_np1;
  	}
  
	void InitObj(){									//Initial state. Not really rendered in Houdini...Just for record
		T a=(T) 0;
    	char str[12];
    	sprintf(str, "%d", 0);
    	std::string frame_name(str);

    	//std::string positions_filename(std::string("particle_x_")+frame_name);
		std::ofstream outfile(std::string("/Users/ninjacat/Desktop/UCLA16FALL/270A/git/Math_270A_Fall_2016_HW3/obj/")+std::string("frame_")+ frame_name + std::string(".obj"));
		
		for(int e=0;e<=N-1;e++){
			a=(T).1/(T)0.7;												//epsilon=0.1
			outfile << "v " << x_n(e) << " " << a << " " << a << " " <<std::endl;
			outfile << "v " << x_n(e) << " " << a << " " << -a << " " <<std::endl;
			outfile << "v " << x_n(e) << " " << -a << " " << a << " " <<std::endl;
			outfile << "v " << x_n(e) << " " << -a << " " << -a << " " <<std::endl;
		}
		
		//e=0;
		outfile << "f " << 1 << " " << 2 << " " << 4 << " " <<std::endl;
		outfile << "f " << 1 << " " << 3 << " " << 4 << " " <<std::endl;
		
		for(int e=1;e<=N-1;e++){
			outfile << "f " << 4*(e-1)+1 << " " << 4*(e-1)+2 << " " << 4*e+2 << " " <<std::endl;
			outfile << "f " << 4*(e-1)+1 << " " << 4*e+1 << " " << 4*e+2 << " " <<std::endl;
			
			outfile << "f " << 4*(e-1)+1 << " " << 4*(e-1)+3 << " " << 4*e+1 << " " <<std::endl;
			outfile << "f " << 4*(e-1)+3 << " " << 4*e+1 << " " << 4*e+3 << " " <<std::endl;
			
			outfile << "f " << 4*(e-1)+3 << " " << 4*e << " " << 4*e+3 << " " <<std::endl;
			outfile << "f " << 4*e+3 << " " << 4*e << " " << 4*(e+1) << " " <<std::endl;
			
			outfile << "f " << 4*(e-1)+2 << " " << 4*e << " " << 4*(e+1) << " " <<std::endl;
			outfile << "f " << 4*(e-1)+2 << " " << 4*e+2 << " " << 4*(e+1) << " " <<std::endl;
		}
		
		//e=N-1;
		outfile << "f " << 4*(N-1)+1 << " " << 4*(N-1)+2 << " " << 4*(N-1)+4 << " " <<std::endl;
		outfile << "f " << 4*(N-1)+1 << " " << 4*(N-1)+3 << " " << 4*(N-1)+4 << " " <<std::endl;
		
		outfile.close();
		
	}
  
  
	void WriteObj(const int number){
		T a=(T) 0;
    	char str[12];
    	sprintf(str, "%d", number);
    	std::string frame_name(str);

    	//std::string positions_filename(std::string("particle_x_")+frame_name);
		std::ofstream outfile(std::string("/Users/ninjacat/Desktop/UCLA16FALL/270A/git/Math_270A_Fall_2016_HW3/obj/")+std::string("frame_")+ frame_name + std::string(".obj"));
		
		//a=(T).1/(T)0.7;															// the left surface is glued to the wall
		a=(T).1/std::sqrt(((T) x_np1(1)-x_np1(0))/( (T) dX));						// only the centerpoint is glued
		outfile << "v " << x_np1(0) << " " << a << " " << a << " " <<std::endl;
		outfile << "v " << x_np1(0) << " " << a << " " << -a << " " <<std::endl;
		outfile << "v " << x_np1(0) << " " << -a << " " << a << " " <<std::endl;
		outfile << "v " << x_np1(0) << " " << -a << " " << -a << " " <<std::endl;
		
		for(int e=1;e<=N-1;e++){
			a=(T).1/std::sqrt(((T) x_np1(e)-x_np1(e-1))/( (T) dX));					//epsilon=0.1
			outfile << "v " << x_np1(e) << " " << a << " " << a << " " <<std::endl;
			outfile << "v " << x_np1(e) << " " << a << " " << -a << " " <<std::endl;
			outfile << "v " << x_np1(e) << " " << -a << " " << a << " " <<std::endl;
			outfile << "v " << x_np1(e) << " " << -a << " " << -a << " " <<std::endl;
		}
		
		//e=0;
		outfile << "f " << 1 << " " << 2 << " " << 4 << " " <<std::endl;
		outfile << "f " << 1 << " " << 3 << " " << 4 << " " <<std::endl;
		
		for(int e=1;e<=N-1;e++){
			outfile << "f " << 4*(e-1)+1 << " " << 4*(e-1)+2 << " " << 4*e+2 << " " <<std::endl;
			outfile << "f " << 4*(e-1)+1 << " " << 4*e+1 << " " << 4*e+2 << " " <<std::endl;
			
			outfile << "f " << 4*(e-1)+1 << " " << 4*(e-1)+3 << " " << 4*e+1 << " " <<std::endl;
			outfile << "f " << 4*(e-1)+3 << " " << 4*e+1 << " " << 4*e+3 << " " <<std::endl;
			
			outfile << "f " << 4*(e-1)+3 << " " << 4*e << " " << 4*e+3 << " " <<std::endl;
			outfile << "f " << 4*e+3 << " " << 4*e << " " << 4*(e+1) << " " <<std::endl;
			
			outfile << "f " << 4*(e-1)+2 << " " << 4*e << " " << 4*(e+1) << " " <<std::endl;
			outfile << "f " << 4*(e-1)+2 << " " << 4*e+2 << " " << 4*(e+1) << " " <<std::endl;
		}
		
		//e=N-1;
		outfile << "f " << 4*(N-1)+1 << " " << 4*(N-1)+2 << " " << 4*(N-1)+4 << " " <<std::endl;
		outfile << "f " << 4*(N-1)+1 << " " << 4*(N-1)+3 << " " << 4*(N-1)+4 << " " <<std::endl;
		
		outfile.close();
	}

  	void Write_State(const int number,std::string& simulation_data_filename){
    	std::ofstream basic_outdata;
    	basic_outdata.open(simulation_data_filename.c_str(), std::ofstream::out | std::ofstream::app);
    	basic_outdata << N << std::endl;
    	basic_outdata.close();

    	char str[12];
    	sprintf(str, "%d", number);
    	std::string frame_name(str);

    	std::string positions_filename(std::string("particle_x_")+frame_name);
    	FILE_IO::Write_Binary(output_directory,positions_filename,x_n);
    	std::string velocities_filename(std::string("particle_v_")+frame_name);
    	FILE_IO::Write_Binary(output_directory,velocities_filename,v_n);
  	}

  	void Read_State(const int number){
    	Read_State(x_n,v_n,N,output_directory);
    	x_np1.resize(N);residual.resize(N);x_hat.resize(N);
  	}

  	static bool Read_State(TVect& x,TVect& v,int& N,std::string output_directory,const int number){
    	std::ifstream basic_indata;
    	int last_frame;
    	std::string simulation_data_filename(output_directory+std::string("/simulation_info.dat"));
    	basic_indata.open(simulation_data_filename.c_str());
    	basic_indata >> last_frame;
    	basic_indata >> N;
    	basic_indata.close();

    	if(number>last_frame || number<0) 
			return false;

    	char str[12];
    	sprintf(str, "%d", number);
    	std::string frame_name(str);
    	x.resize(N);v.resize(N);
    	std::string positions_filename(std::string("particle_x_")+frame_name);
    	FILE_IO::Read_Binary(output_directory,positions_filename,x);
    	std::string velocities_filename(std::string("particle_v_")+frame_name);
    	FILE_IO::Read_Binary(output_directory,velocities_filename,v);

    	return true;
  	}

};
}
#endif
