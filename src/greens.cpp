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



void green_func::set_tstp(int tstp, const green_func_tstp &G){
  assert(tstp>=-1 && tstp<=nt_ && tstp==G.tstp() && G.size1()==size1_ && G.ntau()==ntau_);
  if(tstp == -1){
    memcpy(mat_,G.matptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
  }
  else{
    memcpy(this->retptr(tstp,0),G.retptr(0),sizeof(cplx)*(tstp+1)*element_size_);
    memcpy(this->tvptr(tstp,0),G.tvptr(0),sizeof(cplx)*(ntau_+1)*element_size_);
    memcpy(this->lesptr(0,tstp),G.lesptr(0),sizeof(cplx)*(tstp+1)*element_size_);
  }
}


void green_func::set_tstp_zero(int tstp){
  assert(tstp>=-1 && tstp <=nt_);
  if(tstp==-1){
    memset(matptr(0), 0, sizeof(cplx)*(ntau_+1)*element_size_);
  }
  else{
    memset(retptr(tstp,0), 0, sizeof(cplx)*(tstp+1)*element_size_);
    memset(tvptr(tstp,0), 0, sizeof(cplx)*(ntau_+1)*element_size_);
    memset(lesptr(0,tstp), 0, sizeof(cplx)*(tstp+1)*element_size_);
  }
}


void green_func::set_tstp(int tstp, const green_func &G){
  assert(tstp>=-1 && tstp<=nt_ && tstp<=G.nt() && G.size1()==size1_ && G.ntau()==ntau_);
  if(tstp == -1){
    memcpy(mat_,G.mat_,sizeof(cplx)*(ntau_+1)*element_size_);
  }
  else{
    memcpy(this->retptr(tstp,0),G.retptr(tstp,0),sizeof(cplx)*(tstp+1)*element_size_);
    memcpy(this->tvptr(tstp,0),G.tvptr(tstp,0),sizeof(cplx)*(ntau_+1)*element_size_);
    memcpy(this->lesptr(0,tstp),G.lesptr(0,tstp),sizeof(cplx)*(tstp+1)*element_size_);
  }
}


void green_func::get_tstp(int tstp, green_func_tstp &G) const {
  assert(tstp>=-1 && tstp<=nt_ && tstp==G.tstp() && G.size1()==size1_ && G.ntau()==ntau_);
  if(tstp == -1){
    memcpy(G.matptr(0),mat_,sizeof(cplx)*(ntau_+1)*element_size_);
  }
  else{
    memcpy(G.retptr(0),this->retptr(tstp,0),sizeof(cplx)*(tstp+1)*element_size_);
    memcpy(G.tvptr(0),this->tvptr(tstp,0),sizeof(cplx)*(ntau_+1)*element_size_);
    memcpy(G.lesptr(0),this->lesptr(0,tstp),sizeof(cplx)*(tstp+1)*element_size_);
  }
}

//======================
//MULTIPLICATIONS
//======================

void green_func::smul(int tstp, cplx weight){
  assert(tstp>=-1 && tstp<=nt_);
  if(tstp==-1){
    int len = (ntau_+1)*element_size_;
    for(int i=0;i<len;i++) mat_[i]*=weight;
  }
  else{
    int len = (ntau_+1)*element_size_;
    cplx *tv = this->tvptr(tstp,0);
    for(int i=0;i<len;i++) tv[i]*=weight;

    len = (tstp+1)*element_size_;
    cplx *les = this->lesptr(0,tstp);
    cplx *ret = this->retptr(tstp,0);
    for(int i=0;i<len;i++){
      les[i]*=weight;
      ret[i]*=weight;
    }
  }
}


void green_func::right_multiply(int tstp, cplx *f0, cplx *ft, cplx weight){
  assert(tstp>=-1 && tstp<=nt_);
  int m;
  cplx *x0, *xtmp, *ftmp;
  xtmp = new cplx[element_size_];
  if(tstp == -1){
    // Matsubara component
    x0=mat_;
    for(m=0;m<=ntau_;m++){
      element_mult(size1_,xtmp,x0,f0);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
    }
  }
  else{
    // TV component
    x0=this->tvptr(tstp,0);
    for(m=0;m<=ntau_;m++){
      element_mult(size1_,xtmp,x0,f0);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
    }
    // Ret component
    x0=this->retptr(tstp,0);
    ftmp=ft;
    for(m=0;m<=tstp;m++){
      element_mult(size1_,xtmp,x0,ftmp);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
      ftmp+=element_size_;
    }
    // Les component
    x0=this->lesptr(0,tstp);
    ftmp=ft+tstp*element_size_;
    for(m=0;m<=tstp;m++){
      element_mult(size1_,xtmp,x0,ftmp);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
    }
  }
  delete[] xtmp;
}


void green_func::right_multiply(int tstp, const function &ft, cplx weight){
  assert(tstp>=-1 && tstp<=nt_);
  assert(ft.nt()>=tstp);
  assert(size1_==ft.size1());
  this->right_multiply(tstp,ft.ptr(-1),ft.ptr(0),weight);
}


void green_func::left_multiply(int tstp, cplx *f0, cplx *ft, cplx weight){
  assert(tstp>=-1 && tstp<=nt_);
  int m;
  cplx *x0, *ftmp, *xtmp;
  
  xtmp = new cplx[element_size_];

  if(tstp==-1){
    x0=mat_;
    for(m=0;m<=ntau_;m++){
      element_mult(size1_,xtmp,x0,f0);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
    }
  }
  else{
    // TV component
    x0=this->tvptr(tstp,0);
    ftmp = ft+tstp*element_size_;
    for(m=0;m<=ntau_;m++){
      element_mult(size1_,xtmp,ftmp,x0);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
    }
    // Ret component
    x0=this->retptr(tstp,0);
    for(m=0;m<=tstp;m++){
      element_mult(size1_,xtmp,ftmp,x0);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
    }
    // Les component
    x0=this->lesptr(0,tstp);
    for(m=0;m<=tstp;m++){
      element_mult(size1_,xtmp,ft+m*element_size_,x0);
      element_smul(size1_,xtmp,weight);
      element_set(size1_,x0,xtmp);
      x0+=element_size_;
    } 
  }
  delete[] xtmp;
}


void green_func::left_multiply(int tstp, const function &ft, cplx weight){
  assert(tstp>=-1 && tstp<=nt_);
  assert(ft.nt()>=tstp);
  assert(size1_==ft.size1());
  this->left_multiply(tstp,ft.ptr(-1),ft.ptr(0),weight);
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

double distance_norm2(int tstp, GREEN &G1, GREEN &G2){
  assert(G1.size1()==G2.size1());
  assert(G1.ntau()==G2.ntau());
  assert(tstp<=G1.nt());
  assert(tstp<=G2.nt());

  int size1 = G1.size1(), ntau=G1.ntau(),i;
  double err = 0;
  cplx *tmp = new cplx[size1*size1];
  if(tstp == -1){
    for(i=0;i<=ntau;i++){
      element_set(size1,tmp,G1.matptr(i));
      element_incr(size1,tmp,-1.,G2.matptr(i));
      err += element_norm2(size1,tmp);
    }
  }
  else{
    for(i=0;i<=tstp;i++){
      element_set(size1,tmp,G1.retptr(tstp,i));
      element_incr(size1,tmp,-1.,G2.retptr(tstp,i));
      err += element_norm2(size1,tmp);
    }
    for(i=0;i<=ntau;i++){
      element_set(size1,tmp,G1.tvptr(tstp,i));
      element_incr(size1,tmp,-1.,G2.tvptr(tstp,i));
      err += element_norm2(size1,tmp);
    }
    for(i=0;i<=tstp;i++){
      element_set(size1,tmp,G1.lesptr(i,tstp));
      element_incr(size1,tmp,-1.,G2.lesptr(i,tstp));
      err += element_norm2(size1,tmp);
    }
  }
  delete[] tmp;
  return err;
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

void tti_green_func::set_tstp_zero(int tstp){
  assert(tstp>=-1 && tstp <=nt_);
  if(tstp==-1){
    memset(matptr(0), 0, sizeof(cplx)*(ntau_+1)*element_size_);
  }
  else{
    memset(retptr(tstp), 0, sizeof(cplx)*element_size_);
    memset(tvptr(tstp,0), 0, sizeof(cplx)*(ntau_+1)*element_size_);
    memset(lesptr(-tstp), 0, sizeof(cplx)*element_size_);
  }
}

}//Namespace
#endif
