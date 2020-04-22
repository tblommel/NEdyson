#ifndef PARAM
#define PARAM

#include <fstream>
#include <cstring>

template<typename T> void find_param(char *file, const  char *flag, T &x){
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

#endif
