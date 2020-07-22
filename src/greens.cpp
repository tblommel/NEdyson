#ifndef GREENS_FUNCTIONS_IMPL
#define GREENS_FUNCITONS_IMPL

#include "greens.h"
namespace NEdyson{

//===========================================
//CONSTRUCTORS
//===========================================


green_func::green_func(){
  les_=0;
  ret_=0;
  tv_=0;
  mat_=0;
  ntau_=0;
  nt_=0;
  size1_=0;
  element_size_=0;
  sig_=-1;
}


green_func::~green_func(){
  delete[] les_;
  delete[] ret_;
  delete[] tv_;
  delete[] mat_;
}


green_func::green_func(int nt, int ntau, int size1, int sig){
  assert(size1>0 && nt>=-1 && ntau>=0 && sig*sig==1);
  nt_=nt;
  ntau_=ntau;
  sig_=sig;
  size1_=size1;
  element_size_=size1*size1;
  mat_ = new cplx[(ntau_+1)*element_size_];
  if(nt_ >= 0 ){
    les_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
    ret_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
  } else {
    les_=0;
    ret_=0;
    tv_=0;
  }
}

green_func::green_func(const green_func &g){
  nt_=g.nt_;
  ntau_=g.ntau_;
  sig_=g.sig_;
  size1_=g.size1_;
  element_size_=size1_*size1_;
  
  mat_ = new cplx[(ntau_+1)*element_size_];
  memcpy(mat_,g.mat_,sizeof(cplx)*(ntau_+1)*element_size_);
  if(nt_ >= 0 ){
    les_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
    ret_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
    memcpy(les_,g.les_,sizeof(cplx)*((nt_+1)*(nt_+2))/2*element_size_);
    memcpy(ret_,g.ret_,sizeof(cplx)*((nt_+1)*(nt_+2))/2*element_size_);
    memcpy(tv_,g.tv_,sizeof(cplx)*((nt_+1)*(ntau_+1))*element_size_);
  } else {
    les_=0;
    ret_=0;
    tv_=0;
  }
}

green_func& green_func::operator=(const green_func &g){
  if (this == &g) return *this;
  sig_=g.sig_;
  if(nt_ != g.nt_ || ntau_ != g.ntau_ || size1_ != g.size1_){
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    delete[] mat_;
    nt_ = g.nt_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    element_size_ = size1_*size1_;
    mat_ = new cplx[(ntau_+1)*element_size_];
    memcpy(mat_,g.mat_,sizeof(cplx)*(ntau_+1)*element_size_);
    if(nt_ >= 0){
      les_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
      ret_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
      tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
    } else {
      les_=0;
      ret_=0;
      tv_=0;
    }
  }
  memcpy(les_,g.les_,sizeof(cplx)*((nt_+1)*(nt_+2))/2*element_size_);
  memcpy(ret_,g.ret_,sizeof(cplx)*((nt_+1)*(nt_+2))/2*element_size_);
  memcpy(tv_,g.tv_,sizeof(cplx)*((nt_+1)*(ntau_+1))*element_size_);
  return *this;
}

//======================
//RESIZE
//======================
void green_func::resize(int nt, int ntau, int size1){
  assert(ntau >= 0 && nt >= -1 && size1 > 0);
  delete[] les_;
  delete[] ret_;
  delete[] tv_;
  delete[] mat_;

  nt_=nt;
  ntau_=ntau;
  size1_=size1;
  element_size_=size1*size1;
  mat_ = new cplx[(ntau_+1)*element_size_];
  if(nt_ >= 0 ){
    les_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
    ret_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
  } else {
    les_=0;
    ret_=0;
    tv_=0;
  }
}



//======================
//GET POINTERS
//======================

cplx *green_func::lesptr(int i, int j) const {
  assert(i>=0 && j>= 0 && i<=j && j<=nt_);
  return les_ + ((j*(j+1))/2+i)*element_size_;
}

cplx *green_func::retptr(int i, int j) const {
  assert(i>=0 && j>= 0 && i<=nt_ && j<=i);
  return ret_ + ((i*(i+1))/2+j)*element_size_;
}

cplx *green_func::tvptr(int i, int j) const {
  assert(i>=0 && j>= 0 && i<=nt_ && j<=ntau_);
  return tv_ + (i*(ntau_+1)+j)*element_size_;
}

cplx *green_func::matptr(int i) const {
  assert(i>=0 && i<=ntau_);
  return mat_ + i*element_size_;
}



//======================
//INPUT OUTPUT
//======================
void green_func::print_to_file_mat(std::string file, double dt, double dtau, int precision) const {
  int i, j, l, sg=element_size_,t,tp;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GM.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(i=0;i<=ntau_;i++){
    x=matptr(i);
    for(j=0;j<sg;j++){
      out<<x[j].real()<<" "<<x[j].imag()<<" ";
    }
    out<<std::endl;
  }
  out.close();
}

void green_func::print_to_file_ret(std::string file, double dt, double dtau, int precision) const {
  int i, t, tp, sg = element_size_;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GR.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(t=0;t<=nt_;t++){
    for(tp=0;tp<=t;tp++){
      x=retptr(t,tp);
      for(i=0;i<sg;i++){
        out<<x[i].real()<<" "<<x[i].imag()<<" ";
      }
      out<<std::endl;
    }
  }
  out.close();
}

void green_func::print_to_file_les(std::string file, double dt, double dtau, int precision) const {
  int i, j, l, sg = element_size_;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GL.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(j=0;j<=nt_;j++){
    for(i=0;i<=j;i++){
      x=lesptr(i,j);
      for(l=0;l<sg;l++){
        out<<x[l].real()<<" "<<x[l].imag()<<" ";
      }
      out<<std::endl;
    }
  }
  out.close();
}

void green_func::print_to_file_tv(std::string file, double dt, double dtau, int precision) const {
  int i, j, l, sg = element_size_;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GTV.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(i=0;i<=nt_;i++){
    for(j=0;j<=ntau_;j++){
      x=tvptr(i,j);
      for(l=0;l<sg;l++){
        out<<x[l].real()<<" "<<x[l].imag()<<" ";
      }
      out<<std::endl;
    }
  }
  out.close();
}


void green_func::print_to_file(std::string file, double dt, double dtau, int precision) const {
  print_to_file_mat(file, dt, dtau, precision);
  if(nt_>=0){
    print_to_file_ret(file, dt, dtau, precision);
    print_to_file_tv(file, dt, dtau, precision);
    print_to_file_les(file, dt, dtau, precision);
  } 
}


void green_func::read_from_file_ret(const char *file, double &dt, double &dtau){
  int i,j,l,size1,sg,sig,nt,ntau;
  double real, imag, dt2, dtau2;
  cplx *x=NULL;

  std::ifstream in;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GR.dat";
  in.open(actfile);
  if(!(in>>nt>>ntau>>size1>>dt>>dtau>>sig)){
    std::cerr << "Error in reading greens func from file " << actfile << std::endl;
    abort();
  }

  if(nt!=nt_ || ntau!=ntau_ || size1!=size1_) resize(nt, ntau, size1);
  sig_=sig;
  sg=element_size_;

  delete[] ret_;
  ret_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];

  for(i=0;i<=nt_;i++){
    for(j=0;j<=i;j++){
      x=retptr(i,j);
      for(l=0;l<sg;l++){
        if(!(in>>real>>imag)){
          std::cerr << "Error in reading greens function.  Error in Ret in "<<actfile<<std::endl;
        }
        x[l]=std::complex<double>(real,imag);
      }
    }
  }
  in.close();
}


void green_func::read_from_file(const char *file, double &dt, double &dtau){
  int i,j,l,size1,sg,sig,nt,ntau;
  double real, imag, dt2, dtau2;
  cplx *x=NULL;

  std::ifstream in;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GM.dat";
  in.open(actfile);

  if(!(in>>nt>>ntau>>size1>>dt>>dtau>>sig)){
    std::cerr << "Error in reading greens func from file " << actfile << std::endl;
    abort();
  }

  if(nt!=nt_ || ntau!=ntau_ || size1!=size1_) resize(nt, ntau, size1);
  sig_=sig;
  sg=element_size_;

  for(j=0;j<=ntau_;j++){
    x=matptr(j);
    for(l=0;l<sg;l++){
      if(!(in>>real>>imag)){
        std::cerr << "Error in reading greens function. Error in matsubara in "<<actfile<<std::endl;
      }
      x[l]=std::complex<double>(real,imag);
    }
  }
  in.close();

  if(nt>=0){
    actfile = "";
    actfile += file;
    actfile += "_GR.dat";
    in.open(actfile);
    if(!(in>>nt>>ntau>>size1>>dt2>>dtau2>>sig)){
      std::cerr << "Error in reading greens func from file " << actfile << std::endl;
      abort();
    }
    if(dtau2!=dtau || dt2!=dt || nt!=nt_ || ntau!=ntau_ || size1!=size1_ || sig!=sig_){
      std::cerr << "Ret file parameters are not consistent " << actfile << std::endl;
      abort();
    }
  
    delete[] ret_;
    ret_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];

    for(i=0;i<=nt_;i++){
      for(j=0;j<=i;j++){
        x=retptr(i,j);
        for(l=0;l<sg;l++){
          if(!(in>>real>>imag)){
            std::cerr << "Error in reading greens function.  Error in Ret in "<<actfile<<std::endl;
          }
          x[l]=std::complex<double>(real,imag);
        }
      }
    }
    in.close();

    actfile = "";
    actfile += file;
    actfile += "_GL.dat";
    in.open(actfile);
    if(!(in>>nt>>ntau>>size1>>dt2>>dtau2>>sig)){
      std::cerr << "Error in reading greens func from file " << actfile << std::endl;
      abort();
    }
    if(dtau2!=dtau || dt2!=dt || nt!=nt_ || ntau!=ntau_ || size1!=size1_ || sig!=sig_){
      std::cerr << "Les file parameters are not consistent " << actfile << std::endl;
      abort();
    }

    delete[] les_;
    les_ = new cplx[((nt_+1)*(nt_+2))/2*element_size_];

    for(j=0;j<=nt_;j++){
      for(i=0;i<=j;i++){
        x=lesptr(i,j);
        for(l=0;l<sg;l++){
          if(!(in>>real>>imag)){
            std::cerr << "Error in reading greens function.  Error in Les in "<<actfile<<std::endl;
          }
          x[l]=std::complex<double>(real,imag);
        }
      }
    }
    in.close();




    actfile = "";
    actfile += file;
    actfile += "_GTV.dat";
    in.open(actfile);
    if(!(in>>nt>>ntau>>size1>>dt2>>dtau2>>sig)){
      std::cerr << "Error in reading greens func from file " << actfile << std::endl;
      abort();
    }
    if(dtau2!=dtau || dt2!=dt || nt!=nt_ || ntau!=ntau_ || size1!=size1_ || sig!=sig_){
      std::cerr << "TV file parameters are not consistent " << actfile << std::endl;
      abort();
    }
    delete[] tv_;
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
    for(i=0;i<=nt_;i++){
      for(j=0;j<=ntau_;j++){
        x=tvptr(i,j);
        for(l=0;l<sg;l++){
          if(!(in>>real>>imag)){
            std::cerr << "Error in reading greens function.  Error in Tv in "<<actfile<<std::endl;
          }
          x[l]=std::complex<double>(real,imag);
        }
      }
    }
    in.close();
  }
}


void green_func::print_to_file_mat(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=(ntau_+1)*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GM");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GM",h5::DataSpace(dims));
  dataset.write(mat_);
  File.flush();
}

void green_func::print_to_file_ret(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=((nt_+1)*(nt_+2))/2*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GR");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GR",h5::DataSpace(dims));
  dataset.write(ret_);
  File.flush();
}

void green_func::print_to_file_les(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=((nt_+1)*(nt_+2))/2*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GL");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GL",h5::DataSpace(dims));
  dataset.write(les_);
  File.flush();
}

void green_func::print_to_file_tv(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=((nt_+1)*(ntau_+1))*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GTV");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GTV",h5::DataSpace(dims));
  dataset.write(tv_);
  File.flush();
}

void green_func::print_to_file(h5e::File &File, std::string path) const {
  green_func::print_to_file_mat(File,path);
  green_func::print_to_file_ret(File,path);
  green_func::print_to_file_les(File,path);
  green_func::print_to_file_tv(File,path);
  h5e::dump(File,path+"/nt",nt_);
  h5e::dump(File,path+"/ntau",ntau_);
  h5e::dump(File,path+"/sig",sig_);
  h5e::dump(File,path+"/nao",size1_);
}

void green_func::read_from_file_tv(h5e::File &File,std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GTV");
  size_t len = dataset.getElementCount();
  delete [] tv_;
  tv_ = new cplx[len];
  dataset.read(tv_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void green_func::read_from_file_ret(h5e::File &File,std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GR");
  size_t len = dataset.getElementCount();
  delete [] ret_;
  ret_ = new cplx[len];
  dataset.read(ret_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void green_func::read_from_file_les(h5e::File &File, std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GL");
  size_t len = dataset.getElementCount();
  delete [] les_;
  les_ = new cplx[len];
  dataset.read(les_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void green_func::read_from_file_mat(h5e::File &File,std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GM");
  size_t len = dataset.getElementCount();
  delete [] mat_;
  mat_ = new cplx[len];
  dataset.read(mat_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void green_func::read_from_file(h5e::File &File, std::string path) {
  green_func::read_from_file_mat(File,path);
  green_func::read_from_file_les(File,path);
  green_func::read_from_file_ret(File,path);
  green_func::read_from_file_tv(File,path);
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

//===========================================
//CONSTRUCTORS
//===========================================


tti_green_func::tti_green_func(){
  les_=0;
  ret_=0;
  tv_=0;
  mat_=0;
  ntau_=0;
  nt_=0;
  size1_=0;
  element_size_=0;
  sig_=-1;
}


tti_green_func::~tti_green_func(){
  delete[] les_;
  delete[] ret_;
  delete[] tv_;
  delete[] mat_;
}


tti_green_func::tti_green_func(int nt, int ntau, int size1, int sig){
  assert(size1>0 && nt>=-1 && ntau>=0 && sig*sig==1);
  nt_=nt;
  ntau_=ntau;
  sig_=sig;
  size1_=size1;
  element_size_=size1*size1;
  mat_ = new cplx[(ntau_+1)*element_size_];
  if(nt_ >= 0 ){
    les_ = new cplx[(nt_+1)*element_size_];
    ret_ = new cplx[(nt_+1)*element_size_];
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
  } else {
    les_=0;
    ret_=0;
    tv_=0;
  }
}

tti_green_func::tti_green_func(const tti_green_func &g){
  nt_=g.nt_;
  ntau_=g.ntau_;
  sig_=g.sig_;
  size1_=g.size1_;
  element_size_=size1_*size1_;
  
  mat_ = new cplx[(ntau_+1)*element_size_];
  memcpy(mat_,g.mat_,sizeof(cplx)*(ntau_+1)*element_size_);
  if(nt_ >= 0 ){
    les_ = new cplx[(nt_+1)*element_size_];
    ret_ = new cplx[(nt_+1)*element_size_];
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
    memcpy(les_,g.les_,sizeof(cplx)*(nt_+1)*element_size_);
    memcpy(ret_,g.ret_,sizeof(cplx)*(nt_+1)*element_size_);
    memcpy(tv_,g.tv_,sizeof(cplx)*((nt_+1)*(ntau_+1))*element_size_);
  } else {
    les_=0;
    ret_=0;
    tv_=0;
  }
}

tti_green_func& tti_green_func::operator=(const tti_green_func &g){
  if (this == &g) return *this;
  sig_=g.sig_;
  if(nt_ != g.nt_ || ntau_ != g.ntau_ || size1_ != g.size1_){
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    delete[] mat_;
    nt_ = g.nt_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    element_size_ = size1_*size1_;
    mat_ = new cplx[(ntau_+1)*element_size_];
    memcpy(mat_,g.mat_,sizeof(cplx)*(ntau_+1)*element_size_);
    if(nt_ >= 0){
      les_ = new cplx[(nt_+1)*element_size_];
      ret_ = new cplx[(nt_+1)*element_size_];
      tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
    } else {
      les_=0;
      ret_=0;
      tv_=0;
    }
  }
  memcpy(les_,g.les_,sizeof(cplx)*(nt_+1)*element_size_);
  memcpy(ret_,g.ret_,sizeof(cplx)*(nt_+1)*element_size_);
  memcpy(tv_,g.tv_,sizeof(cplx)*((nt_+1)*(ntau_+1))*element_size_);
  return *this;
}
//======================
//RESIZE
//======================
void tti_green_func::resize(int nt, int ntau, int size1){
  assert(ntau >= 0 && nt >= -1 && size1 > 0);
  delete[] les_;
  delete[] ret_;
  delete[] tv_;
  delete[] mat_;

  nt_=nt;
  ntau_=ntau;
  size1_=size1;
  element_size_=size1*size1;
  mat_ = new cplx[(ntau_+1)*element_size_];
  if(nt_ >= 0 ){
    les_ = new cplx[(nt_+1)*element_size_];
    ret_ = new cplx[(nt_+1)*element_size_];
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
  } else {
    les_=0;
    ret_=0;
    tv_=0;
  }
}



//======================
//GET POINTERS
//======================

cplx *tti_green_func::lesptr(int i) const {
  assert(i<=0 && i>=-nt_);
  return les_ - i*element_size_;
}

cplx *tti_green_func::retptr(int i) const {
  assert(i>=0 && i<=nt_);
  return ret_ + i*element_size_;
}

cplx *tti_green_func::tvptr(int i, int j) const {
  assert(i>=0 && j>= 0 && i<=nt_ && j<=ntau_);
  return tv_ + (i*(ntau_+1)+j)*element_size_;
}

cplx *tti_green_func::matptr(int i) const {
  assert(i>=0 && i<=ntau_);
  return mat_ + i*element_size_;
}



//======================
//INPUT OUTPUT
//======================
void tti_green_func::print_to_file_mat(std::string file, double dt, double dtau, int precision) const {
  int i, j, l, sg=element_size_,t,tp;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GM.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(i=0;i<=ntau_;i++){
    x=matptr(i);
    for(j=0;j<sg;j++){
      out<<x[j].real()<<" "<<x[j].imag()<<" ";
    }
    out<<std::endl;
  }
  out.close();
}

void tti_green_func::print_to_file_ret(std::string file, double dt, double dtau, int precision) const {
  int i, t, tp, sg = element_size_;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GR.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(t=0;t<=nt_;t++){
    x=retptr(t);
    for(i=0;i<sg;i++){
      out<<x[i].real()<<" "<<x[i].imag()<<" ";
    }
    out<<std::endl;
  }
  out.close();
}

void tti_green_func::print_to_file_les(std::string file, double dt, double dtau, int precision) const {
  int i, j, l, sg = element_size_;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GL.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(j=0;j>=-nt_;j--){
    x=lesptr(j);
    for(l=0;l<sg;l++){
      out<<x[l].real()<<" "<<x[l].imag()<<" ";
    }
    out<<std::endl;
  }
  out.close();
}

void tti_green_func::print_to_file_tv(std::string file, double dt, double dtau, int precision) const {
  int i, j, l, sg = element_size_;
  cplx *x = NULL;
  std::ofstream out;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GTV.dat";
  out.open(actfile);
  out.precision(precision);
  out << nt_ << " " << ntau_ << " " << size1_ << " " << dt <<" "<<dtau<<" " << sig_ << std::endl;
  for(i=0;i<=nt_;i++){
    for(j=0;j<=ntau_;j++){
      x=tvptr(i,j);
      for(l=0;l<sg;l++){
        out<<x[l].real()<<" "<<x[l].imag()<<" ";
      }
      out<<std::endl;
    }
  }
  out.close();
}


void tti_green_func::print_to_file(std::string file, double dt, double dtau, int precision) const {
  print_to_file_mat(file, dt, dtau, precision);
  if(nt_>=0){
    print_to_file_ret(file, dt, dtau, precision);
    print_to_file_tv(file, dt, dtau, precision);
    print_to_file_les(file, dt, dtau, precision);
  } 
}



void tti_green_func::read_from_file_ret(const char *file, double &dt, double &dtau){
  int i,j,l,size1,sg,sig,nt,ntau;
  double real, imag, dt2, dtau2;
  cplx *x=NULL;

  std::ifstream in;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GR.dat";
  in.open(actfile);
  if(!(in>>nt>>ntau>>size1>>dt>>dtau>>sig)){
    std::cerr << "Error in reading greens func from file " << actfile << std::endl;
    abort();
  }

  if(nt!=nt_ || ntau!=ntau_ || size1!=size1_) resize(nt, ntau, size1);
  sig_=sig;
  sg=element_size_;

  delete[] ret_;
  ret_ = new cplx[(nt_+1)*element_size_];

  for(i=0;i<=nt_;i++){
    x=retptr(i);
    for(l=0;l<sg;l++){
      if(!(in>>real>>imag)){
        std::cerr << "Error in reading greens function.  Error in Ret in "<<actfile<<std::endl;
      }
      x[l]=std::complex<double>(real,imag);
    }
  }
  in.close();
}


void tti_green_func::read_from_file(const char *file, double &dt, double &dtau){
  int i,j,l,size1,sg,sig,nt,ntau;
  double real, imag, dt2, dtau2;
  cplx *x=NULL;

  std::ifstream in;
  std::string actfile;
  actfile = "";
  actfile += file;
  actfile += "_GM.dat";
  in.open(actfile);

  if(!(in>>nt>>ntau>>size1>>dt>>dtau>>sig)){
    std::cerr << "Error in reading greens func from file " << actfile << std::endl;
    abort();
  }

  if(nt!=nt_ || ntau!=ntau_ || size1!=size1_) resize(nt, ntau, size1);
  sig_=sig;
  sg=element_size_;

  for(j=0;j<=ntau_;j++){
    x=matptr(j);
    for(l=0;l<sg;l++){
      if(!(in>>real>>imag)){
        std::cerr << "Error in reading greens function. Error in matsubara in "<<actfile<<std::endl;
      }
      x[l]=std::complex<double>(real,imag);
    }
  }
  in.close();

  if(nt>=0){
    actfile = "";
    actfile += file;
    actfile += "_GR.dat";
    in.open(actfile);
    if(!(in>>nt>>ntau>>size1>>dt2>>dtau2>>sig)){
      std::cerr << "Error in reading greens func from file " << actfile << std::endl;
      abort();
    }
    if(dtau2!=dtau || dt2!=dt || nt!=nt_ || ntau!=ntau_ || size1!=size1_ || sig!=sig_){
      std::cerr << "Ret file parameters are not consistent " << actfile << std::endl;
      abort();
    }
  
    delete[] ret_;
    ret_ = new cplx[(nt_+1)*element_size_];

    for(i=0;i<=nt_;i++){
      x=retptr(i);
      for(l=0;l<sg;l++){
        if(!(in>>real>>imag)){
          std::cerr << "Error in reading greens function.  Error in Ret in "<<actfile<<std::endl;
        }
        x[l]=std::complex<double>(real,imag);
      }
    }
    in.close();

    actfile = "";
    actfile += file;
    actfile += "_GL.dat";
    in.open(actfile);
    if(!(in>>nt>>ntau>>size1>>dt2>>dtau2>>sig)){
      std::cerr << "Error in reading greens func from file " << actfile << std::endl;
      abort();
    }
    if(dtau2!=dtau || dt2!=dt || nt!=nt_ || ntau!=ntau_ || size1!=size1_ || sig!=sig_){
      std::cerr << "Les file parameters are not consistent " << actfile << std::endl;
      abort();
    }

    delete[] les_;
    les_ = new cplx[(nt_+1)*element_size_];

    for(j=0;j>=-nt_;j--){
      x=lesptr(j);
      for(l=0;l<sg;l++){
        if(!(in>>real>>imag)){
          std::cerr << "Error in reading greens function.  Error in Les in "<<actfile<<std::endl;
        }
        x[l]=std::complex<double>(real,imag);
      }
    }
    in.close();




    actfile = "";
    actfile += file;
    actfile += "_GTV.dat";
    in.open(actfile);
    if(!(in>>nt>>ntau>>size1>>dt2>>dtau2>>sig)){
      std::cerr << "Error in reading greens func from file " << actfile << std::endl;
      abort();
    }
    if(dtau2!=dtau || dt2!=dt || nt!=nt_ || ntau!=ntau_ || size1!=size1_ || sig!=sig_){
      std::cerr << "TV file parameters are not consistent " << actfile << std::endl;
      abort();
    }
    delete[] tv_;
    tv_ = new cplx[(nt_+1)*(ntau_+1)*element_size_];
    for(i=0;i<=nt_;i++){
      for(j=0;j<=ntau_;j++){
        x=tvptr(i,j);
        for(l=0;l<sg;l++){
          if(!(in>>real>>imag)){
            std::cerr << "Error in reading greens function.  Error in Tv in "<<actfile<<std::endl;
          }
          x[l]=std::complex<double>(real,imag);
        }
      }
    }
    in.close();
  }
}


void tti_green_func::print_to_file_mat(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=(ntau_+1)*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GM");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GM",h5::DataSpace(dims));
  dataset.write(mat_);
  File.flush();
}

void tti_green_func::print_to_file_ret(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=(nt_+1)*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GR");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GR",h5::DataSpace(dims));
  dataset.write(ret_);
  File.flush();
}

void tti_green_func::print_to_file_les(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=(nt_+1)*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GL");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GL",h5::DataSpace(dims));
  dataset.write(les_);
  File.flush();
}

void tti_green_func::print_to_file_tv(h5e::File &File, std::string path) const {
  std::vector<size_t> dims(1);
  dims[0]=((nt_+1)*(ntau_+1))*element_size_;
  h5e::detail::createGroupsToDataSet(File, path+"/GTV");
  h5::DataSet dataset = File.createDataSet<cplx>(path+"/GTV",h5::DataSpace(dims));
  dataset.write(tv_);
  File.flush();
}

void tti_green_func::print_to_file(h5e::File &File, std::string path) const {
  tti_green_func::print_to_file_mat(File,path);
  tti_green_func::print_to_file_ret(File,path);
  tti_green_func::print_to_file_les(File,path);
  tti_green_func::print_to_file_tv(File,path);
  h5e::dump(File,path+"/nt",nt_);
  h5e::dump(File,path+"/ntau",ntau_);
  h5e::dump(File,path+"/sig",sig_);
  h5e::dump(File,path+"/nao",size1_);
}

void tti_green_func::read_from_file_tv(h5e::File &File,std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GTV");
  size_t len = dataset.getElementCount();
  delete [] tv_;
  tv_ = new cplx[len];
  dataset.read(tv_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void tti_green_func::read_from_file_ret(h5e::File &File,std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GR");
  size_t len = dataset.getElementCount();
  delete [] ret_;
  ret_ = new cplx[len];
  dataset.read(ret_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void tti_green_func::read_from_file_les(h5e::File &File, std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GL");
  size_t len = dataset.getElementCount();
  delete [] les_;
  les_ = new cplx[len];
  dataset.read(les_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void tti_green_func::read_from_file_mat(h5e::File &File,std::string path) {
  h5::DataSet dataset = File.getDataSet(path+"/GM");
  size_t len = dataset.getElementCount();
  delete [] mat_;
  mat_ = new cplx[len];
  dataset.read(mat_);
  File.flush();
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

void tti_green_func::read_from_file(h5e::File &File, std::string path) {
  tti_green_func::read_from_file_mat(File,path);
  tti_green_func::read_from_file_les(File,path);
  tti_green_func::read_from_file_ret(File,path);
  tti_green_func::read_from_file_tv(File,path);
  nt_ = h5e::load<int>(File, path+"/nt");
  sig_ = h5e::load<int>(File,path+"/sig");
  ntau_ = h5e::load<int>(File,path+"/ntau");
  size1_ = h5e::load<int>(File,path+"/nao");
  element_size_=size1_*size1_;
}

}//Namespace
#endif
