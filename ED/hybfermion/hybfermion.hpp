#pragma once
#include<Eigen/Dense>
#include<vector>
#include<fstream>

class fermion_op{
public:
  typedef Eigen::Matrix< double , Eigen::Dynamic, Eigen::Dynamic > matrix_t;
  fermion_op(int orbitals):
    orbitals_(orbitals),
    dim_(1<<orbitals_){
    create_creation_ops();
    validate();
    print_to_file();
  }
  int orbital_mask(int i){ return 1<<i;}
  int orbitals() const{return orbitals_;}
  int dim() const{return dim_;}
  int cdag(int orbital, int targ, int src) const{return cdag_ops[orbital](src, targ); }
  int c(int orbital, int targ, int src) const{ return c_ops[orbital](src, targ); }
  void print_to_file(){
    for(int i=0;i<orbitals_;++i){
      std::stringstream cname; cname<<"c_"<<i<<".dat";
      std::stringstream cdagname; cdagname<<"cdag_"<<i<<".dat";
      std::ofstream c_file(cname.str().c_str()); c_file<<c_ops[i];
      std::ofstream cdag_file(cdagname.str().c_str()); cdag_file<<cdag_ops[i];
    }
  }
private:
  void create_creation_ops(){
    for(int i=0;i<orbitals_;++i){
      matrix_t c_i=matrix_t::Zero(dim_,dim_);
      //annihilator will change source state to target state.
      for(int source_state=0;source_state<dim_;++source_state){
        if(source_state&orbital_mask(i)){
          int target_state=source_state^orbital_mask(i);
          int permutation_sign=1;
          for(int k=i+1;k<orbitals_;++k) if(source_state & orbital_mask(k)) permutation_sign*=-1;
          c_i(source_state, target_state)=permutation_sign;
        }
      }
      c_ops.push_back(c_i);
      cdag_ops.push_back(c_i.transpose());
    }
  }
  ///check all anticommutation relations
  void validate() const{
    for(int i=0;i<orbitals_;++i){
      for(int j=0;j<orbitals_;++j){
        if(anti_commutator(c_ops[i], c_ops[j]).maxCoeff()>1.e-5) throw std::runtime_error("fermionic operators do not behave as expected: c ops");
        if(anti_commutator(cdag_ops[i], cdag_ops[j]).maxCoeff()>1.e-5) throw std::runtime_error("fermionic operators do not behave as expected: cdag ops");
        if((anti_commutator(cdag_ops[i], c_ops[j])-((i==j)?matrix_t(matrix_t::Identity(dim_,dim_)):matrix_t(matrix_t::Zero(dim_,dim_)))).maxCoeff()>1.e-5){
          std::cout<<"cdag: "<<cdag_ops[i]<<std::endl;
          std::cout<<"c    : "<<c_ops[i]<<std::endl;
          std::cout<<i<<" "<<j<<std::endl<< anti_commutator(cdag_ops[i], c_ops[j])<<std::endl;
              throw std::runtime_error("fermionic operators do not behave as expected.");
        }
      }
    }
  }
  ///the anticommutation relations of fermion operators
  matrix_t anti_commutator(const matrix_t &i, const matrix_t &j) const{return i*j+j*i;}
  ///creation operators: one per (spin-)orbital
  std::vector<matrix_t> c_ops;
  ///annihilation operators: one per (spin-)orbital
  std::vector<matrix_t> cdag_ops;
  ///number of orbitals
  const int orbitals_;
  ///dimension of fock space
  const int dim_;
};
