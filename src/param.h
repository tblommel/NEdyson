#ifndef PARAM
#define PARAM

#include <fstream>
#include <cstring>
#include <Eigen/Eigen>

void find_param(char *file, const  char *flag, int &x){
	std::string actfile;
	actfile = "";
	actfile += file;
	actfile += "input.inp";
  std::ifstream in(actfile);
  std::string str;
  std::string sflag(flag);
  int found=0;
  while(std::getline(in,str)){
    if(str.find(sflag)!=std::string::npos){
      int p1=str.find("=");
      std::stringstream(str.substr(p1+1,str.length()))>>x;
      found=1;
      in.close();
    }
  }
  if(found==0) throw std::invalid_argument(std::string("could not find ")+ std::string(flag));
}

void find_param(char *file, const  char *flag, double &x){
	std::string actfile;
	actfile = "";
	actfile += file;
	actfile += "input.inp";
  std::ifstream in(actfile);
  std::string str;
  std::string sflag(flag);
  int found=0;
  while(std::getline(in,str)){
    if(str.find(sflag)!=std::string::npos){
      int p1=str.find("=");
      std::stringstream(str.substr(p1+1,str.length()))>>x;
      found=1;
      in.close();
    }
  }
  if(found==0) throw std::invalid_argument(std::string("could not find ")+ std::string(flag));
}

void find_param(char *file, const  char *flag, Eigen::VectorXd &x){
	std::string actfile;
	actfile = "";
	actfile += file;
	actfile += "input.inp";
  std::ifstream in(actfile);
  std::string str;
  std::string sflag(flag);
  int found=0;
  while(std::getline(in,str)){
    if(str.find(sflag)!=std::string::npos){
      int p1=str.find("=");
      std::istringstream iss(str.substr(p1+1,str.length()));
      int i=0;
      while(std::getline(iss, str,' ')){
        std::istringstream val(str);
        val>>x(i);
        i++;
      }
      found=1;
      in.close();
    }
  }
  if(found==0) throw std::invalid_argument(std::string("could not find ")+ std::string(flag));
}

#endif
