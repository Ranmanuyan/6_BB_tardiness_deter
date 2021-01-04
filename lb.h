/**
 *   \file  lb.h
 *   \brief  This file group the classes/functions use to compute the Permutation 
 *   Flow shop Lower Bound
 *
 *  \author  Bertrand Le Cun (blec), Bertrand.Lecun@prism.uvsq.fr
 *
 *  \internal
 *    Created:  12.12.2007
 *   Revision:  $Id: doxygen.templates.example,v 1.4 2007/08/02 14:35:24 mehner Exp $
 *   Compiler:  gcc/g++
 *    Company:  PRiSM Laboratory University of Versailles-Saint-Quentin
 *  Copyright:  Copyright (c) 2007, Bertrand Le Cun
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 * =====================================================================================
 */


#include"pfs.h"

/** Think about  : the OneMachine and TwoMachine objects should be stored on the Instance.
 * then the function LB1 and LB2 take the pfi the per and the Scheduled
 * and compute the bound for this node.
 * At this time the construction of these two objects is made at each evaluation... 
 */

/**
 *  \class Jd 
 *  \brief class Job description that stores a value t and an identifier id.
 */
struct Jd {
   /// the value
   int t;
   /// the identifier
   int id;
   /* Pack method
    */
   void Pack(Bob::Serialize &bs) const {
         Bob::Pack(bs,&id);
         Bob::Pack(bs,&t);
   } 
   /* UnPack method
    */
   void UnPack(Bob::DeSerialize &ds) {
         Bob::UnPack(ds,&id);
         Bob::UnPack(ds,&t);
   } 
};

inline std::ostream &operator<<(std::ostream &os,const Jd &j) {
  os <<"("<<j.t<<","<<j.id<<")";
  return os;
}

/**
 *  \brief compare two Jd objects (according to the t value)
 *  \param  j1 the fist Jd object
 *  \param  j2 the second Jd object
 *  \return a boolean 
 */
bool Jd_cmpC(const Jd &j1,const Jd &j2) {
   return j1.t<j2.t;
}

/**
 *  \brief compare two Jd objects (according to the t value)
 *  \param  j1 the fist Jd object
 *  \param  j2 the second Jd object
 *  \return a boolean 
 */
bool Jd_cmpD(const Jd &j1,const Jd &j2) {
   return j1.t>j2.t;
}


/** type Jsd is an array of Jd
 */
typedef Bob::pvector<Jd> Jsd;

/** type Jssd is an array of Jsd
 * then Jssd is a two dimensional array of Jd
 */
typedef Bob::pvector< Jsd > Jssd;

/** 
 * \brief struct J_Order store a job permutation and the job lag associated to a 
 * couple of machine.
 */
// struct J_Order {
//    Jsd lag;
//    Bob::Permutation pm;

//    J_Order() : lag(),pm() {}
//    J_Order(int n) : lag(n),pm(n) {}
//    void Pack(Bob::Serialize &bs) const {
//          Bob::Pack(bs,&lag);
//          Bob::Pack(bs,&pm);
//    } 
//    void UnPack(Bob::DeSerialize &ds) {
//          Bob::UnPack(ds,&lag);
//          Bob::UnPack(ds,&pm);
//    } 

// };

// typedef Bob::pvector< Bob::pvector< J_Order > > l_Date;


// inline std::ostream &operator<<(std::ostream &os,const J_Order &jo) {
//   jo.pm.Prt(os);
//   return os;
// }




void DispJdTab(Jsd &v,int nb) {
      for (int j=0;j<nb;j++ ) { 
            std::cout<< "("<<v[j].id<<","<<v[j].t<<")";
      }
      std::cout<<std::endl;
}

void ConsRelease(Jssd &r,const PFSInstance *pfi) {

   for (int i=0;i<pfi->nbm;i++ ) {
      r[i].resize(pfi->nbj);
      for (int j=0;j<pfi->nbj;j++ ) { r[i][j].t=0;r[i][j].id=j; }
   }
   for (int i=1;i<pfi->nbm;i++ ) {
      for (int j=0;j<pfi->nbj;j++ ) { 
         r[i][j].t=r[i-1][j].t+pfi->d[j][i-1];
      }
      //std::sort(r[i].begin(),r[i].end(),Jd_cmpC);
   }
}



void ConsDelivery(Jssd &q,const PFSInstance *pfi) {

   for (int i=0;i<pfi->nbm;i++ ) {
      q[i].resize(pfi->nbj);
      for (int j=0;j<pfi->nbj;j++ ) { q[i][j].t=0;q[i][j].id=j; }
   }
   for (int i=pfi->nbm-2;i>-1;i-- ) {
      for (int j=0;j<pfi->nbj;j++ ) { 
          q[i][j].t=q[i+1][j].t+pfi->d[j][i+1];
      }
      //std::sort(q[i].begin(),q[i].end(),Jd_cmpC);
   }
}

/*int getMin(Jsd &r) {
   int id=0;
   for (int i=0;i<r.size();i++) {
        if ( r[i].t<r[id].t ) id=i;
   }
   return id;
}*/

int getMinFree(Jsd &r,const Bob::Permutation &per) {
   for (unsigned int i=0;i<r.size();i++) {
        if ( per.isfree(r[i].id) ) return i;
   }
   return -1;
}

int getMinFreeNS(Jsd &r,const Bob::Permutation &per) {
   int min=-1;
   for (unsigned int i=0;i<r.size();i++) {
        if ( per.isfree(r[i].id) && (min==-1 || r[min].t>r[i].t )) { min=i; }
   }
   return min;
}

int getValMinFreeNS(Jsd &r,const Bob::Permutation &per) {
   int min=getMinFreeNS(r,per);
   if ( min==-1 ) { return 0; }
   return r[min].t;
}

int getpk(int m,const PFSInstance *pfi,const Bob::Permutation &per) {
   int pm=0;
   for (int j=0;j<pfi->nbj;j++ ) { 
      if ( per.isfree(j) ) {
         pm+=pfi->d[j][m];
      }
   }
   return pm;
}

// *
// *
// *
// *  The methods that used to calculate the LB of the total tardiness 
// *  6 June 2020
// *  @ Lei Liu
// *
// *
// *
// /

// partial schedule
// the completion time of each job on each machine.
// finishTime[i][j] i is machine index 
int** calculateC(vector<int> PS, const PFSInstance *pfi){
   
   int**  finishTime = new int* [pfi->nbm];
   for(int i=0;i<pfi->nbm;i++)
      finishTime[i] = new int[PS.size()];

   for(int i=0;i<pfi->nbm;i++){
      for(int j=0; j< PS.size();j++){
         if(j==0 && i == 0){
            finishTime[i][j] = pfi->d[PS.front()][0];
         }else if((j==0) && (i!=0)) {
            finishTime[i][j] = finishTime[i-1][j]+pfi->d[PS.front()][i];
         }else if(j != 0 && i == 0){
            finishTime[i][j] = finishTime[i][j-1]+pfi->d[PS[j]][i];
         }else{
            finishTime[i][j] = max( finishTime[i-1][j],finishTime[i][j-1] )+pfi->d[PS[j]][i];
         }
      }
   }
   // for(int j=0;j<this->m;j++)
   //    delete[] finishTime[j];
   // delete[] finishTime;
   PS.clear();
   return finishTime;
}

// partial schedule
// the completion time on each machine of the last job
int* computeL(vector<int> PS, const PFSInstance *pfi){

   int* lastJonC = new int [pfi->nbm];
   int** finishTime = calculateC(PS, pfi);

   if(PS.size() == 0){
      for(int i=0;i<pfi->nbm;i++){
         lastJonC[i] = 0;
      }
   }else if(PS.size() == 1){
      lastJonC[0] = pfi->d[PS.front()][0];
      for(int i=1; i < pfi->nbm; i++){
         lastJonC[i] = lastJonC[i-1]+pfi->d[PS.front()][i];
      }
   }
   else{
      for(int i=0; i < pfi->nbm; i++){
         lastJonC[i] =  finishTime[i][PS.size()-1];
      }
    }

   for(int i=0;i<pfi->nbm;i++)
      delete[] finishTime[i];
   delete[] finishTime;



    return lastJonC;
}

// get the  job index of the unassigned jobs.

vector<int> calNPS(const PFSInstance *pfi, const Bob::Permutation &per){

   vector<int>NPS;
   for(int i=0;i<pfi->nbj;i++) {
      if(per.isfree(i))
        NPS.push_back(i);
    }

    return NPS;
}


// smallest P
// from the unassigned jobs set， sort processing time，根据输入的t, 看取前多少个进行相加。


int smallestP(int k, int t, vector<int> PS,  vector<int> NPS, const PFSInstance *pfi){

   // new vector to store the pt of the NPS on machine k

   vector<int>pNPS;
   for(int i=0; i < NPS.size();i++ ){
      pNPS.push_back(pfi->d[NPS[i]][k]);
   }
   sort(pNPS.begin(),pNPS.end());
   int sumP = 0;
   for(int i=0; i<= t; i++){

      sumP += pNPS[i];

   }
   // clear or delete?
   NPS.clear();
   pNPS.clear();


   return sumP;
}


// smallest value of q(j,r,k-1), r<=k-1, k=1,..,m-1.
// 先求q(j,r,k-1)  得到一个j based 一维数组 qArray = qVector

int smallQ(int r, int k, vector<int>NPS, const PFSInstance *pfi){

   int* qArray = new int[pfi->nbj];
   vector<int>qVecotr;


   for(int j=0; j < NPS.size(); j++){
      for(int u = r; u < k; u++){
         qArray[j] += pfi->d[NPS[j]][u];
      }
      qVecotr.push_back(qArray[j]);
   }

   sort(qVecotr.begin(), qVecotr.end());

   int smallQ = qVecotr.front();

   qVecotr.clear();
   delete []qArray;

   return smallQ;

}


// earliest start time of the t-th job in NPS on machine k.
// two-demension array : t = (0, n-PS.size()-1) , k = (0,m)

// Rk -   one demension
// Rst/EST-  two demension

int** calculateEST(vector<int> PS, const PFSInstance *pfi, const Bob::Permutation &per){
   const vector<int> NPS  = calNPS(pfi, per);
   int**  Est = new int* [pfi->nbj - PS.size()];
   for (int i = 0; i < pfi->nbj - PS.size(); i++) {
      Est[i] = new int[pfi->nbm];
   }
   // Machine 0
   int*  lastJonC = new int [pfi->nbm];
   lastJonC = computeL(PS, pfi);
   Est[0][0] = lastJonC[0];
   int* Rk = new int[pfi->nbm];
   int**  Rst = new int* [pfi->nbj - PS.size()];
   for (int i = 0; i < pfi->nbj - PS.size(); i++) {
      Rst[i] = new int[pfi->nbm];
   }
   Rk[0] = Est[0][0];

   for(int t= 0; t< pfi->nbj - PS.size(); t++){
      Rst[t][0] =  Rk[0] + smallestP(0, t, PS, NPS, pfi);
   }
   for(int t= 1; t< pfi->nbj - PS.size(); t++){
      Est[t][0] = Rst[t-1][0];
   }
   // From  Machine 1 to Machine m-1
   for(int k = 1; k < pfi->nbm; k++){
      int maxr=0;       //EST[0][k]
      for(int r=0; r<k;r++){
         if((Est[0][r]+smallQ(r,k,NPS, pfi)) > maxr){
            maxr = Est[0][r]+smallQ(r,k,NPS, pfi);
         }
      }
      Est[0][k] = max(computeL(PS, pfi)[k], maxr);
      Rk[k] = Est[0][k];
   }
   //Rst[t][k]   Est[t][k]
   for(int k = 1; k < pfi->nbm; k++){
      for(int t= 0; t< pfi->nbj - PS.size(); t++){
         Rst[t][k] = Rk[k] + smallestP(k,t,PS,NPS, pfi);
         //Est[t+1][k] = max(Rst[t][k], Rst[t+1][k-1]);
      }
      for(int t= 1; t< pfi->nbj - PS.size(); t++){
        Est[t][k] = max(Rst[t-1][k], Rst[t][k-1]);
      }
   }
   for(int i = 0; i < pfi->nbj - PS.size(); i++){
      delete[] Rst[i];
   }
   delete[] Rst;
   delete[] lastJonC;
   delete[] Rk;


   return Est;
}


// h function
// chongxie q(jkm)
int** calculateH(const PFSInstance *pfi){

   int**  hValue = new int* [pfi->nbm];
   for (int i = 0; i < pfi->nbm; i++) {
        hValue[i] = new int[pfi->nbj];
    }

   for(int j = 0; j < pfi->nbj; j++){
      for(int k = 0; k < pfi->nbm; k++){
         int q=0;
         for(int r=k; r<pfi->nbm; r++){
            q += pfi->d[j][r];
         }
         hValue[k][j] = pfi->dueDate[j] - q;
      }

   }
   return hValue;
}

//if element in the vector
bool is_element_in_vector(vector<int> v,int element){
   vector<int>::iterator it;
   it=find(v.begin(),v.end(),element);
   if (it!=v.end()){
      return true;
   }
   else{
      return false;
   }
}


//  the t-th smallest value of hValue of the NPS
int smallH(int k, int t, vector<int> PS, const PFSInstance *pfi){

   int **hValue = calculateH(pfi);

   vector<int>qVecotr;

   for(int j=0; j < pfi->nbj; j++){
      if(! is_element_in_vector(PS,j)){
         qVecotr.push_back(hValue[k][j]);
      }

   }

   sort(qVecotr.begin(), qVecotr.end());

   // the t-th smallest element
   // Access element at index 3
   int &element = qVecotr[t];

   qVecotr.clear();

   // free(hValue);

   for(int i = 0; i < pfi->nbm; i++){
      delete[] hValue[i];
   }
   delete[] hValue;

   return element;


}

// use the finish time to calculate total tardiness for jobs in the partial schedule r
int calTardiness(vector<int> PS, const PFSInstance *pfi){
   int**  finishTime = calculateC(PS, pfi);

   //job j  finishTime[m-1][j]
   // job j tardiness  -due[j]
   // sum

   int totalTardiness = 0;

   for(int j=0; j<PS.size(); j++){
      
      totalTardiness += max(finishTime[pfi->nbm -1][j] - pfi->dueDate[PS[j]], 0);

   }

   // run this delete for 9 times, PS.size()= 20, throw error

   for(int i = 0; i < pfi->nbm; i++){
      delete[] finishTime[i];
   }
   delete[] finishTime;

   // free(finishTime);


   return totalTardiness;
}







/**
 *  \class One Machine is a simple lower bound for the one machine problem.
 *  \brief
 *   
 *  \par 
 */
class OneMachine {
   public:
   Jssd r,d;
   Bob::pvector<int> pk;
   int nbj;
   OneMachine() : r(),d(),pk(),nbj(0) { }
   OneMachine(const PFSInstance *pfi) : r(pfi->nbm),d(pfi->nbm),pk(pfi->nbm) {
      ConsRelease(r,pfi);
      ConsDelivery(d,pfi);
      for (int i=0;i<pfi->nbm;i++ ) {
         pk[i]=0;
      }
      for (int i=0;i<pfi->nbm;i++ ) {
         for (int j=0;j<pfi->nbj;j++ ) { 
            pk[i]+=pfi->d[j][i];
         }
      }
   }

   int LB2(const PFSInstance *pfi,const Bob::Permutation &per,Scheduled &md){
      vector<int> j2i =  per.get_j2i();
      vector<int> PS = vector<int>(j2i.begin(), j2i.begin()+per.size()-per.nbFree());
      int** Est = calculateEST(PS, pfi, per);
      int** Tst = new int* [pfi->nbj - PS.size()];
      for (int i = 0; i < (pfi->nbj - PS.size()); i++) {
           Tst[i] = new int[pfi->nbm];
       }
      for(int t=0; t< (pfi->nbj - PS.size()); t++){
         for(int k=0; k<pfi->nbm;k++){
            Tst[t][k] = max(Est[t][k]- smallH( k, t, PS, pfi) , 0);
         }
      }
      int* g = new int [pfi->nbm]();
      
      for(int k=0; k<pfi->nbm;k++){
         for(int t=0; t < (pfi->nbj - PS.size()); t++){
            g[k] += Tst[t][k];
         }
         // std::cout << "g[" << k << "]= " <<g[k] <<endl; 
         // if(g[k] > maxgk) maxgk = g[k];
      }
      int maxgk = 0;
      for(int k=0; k<pfi->nbm;k++){
         if(g[k] > maxgk){
            maxgk = g[k];
         }
      }



      int LB = 0;
      LB = calTardiness(PS, pfi) + maxgk;
      for(int i = 0; i < pfi->nbj - PS.size(); i++){
         delete[] Est[i];
      }
      delete[] Est;
      for(int i = 0; i < pfi->nbj - PS.size(); i++){
         delete[] Tst[i];
      }
      delete[] Tst;

      delete []g;

      j2i.clear();
      PS.clear();

      return LB;
   }

   void Prt(std::ostream &os) {
      os<<" release :"<<std::endl;
      for (unsigned int i=0;i<r.size();i++ ) {
         for (int j=0;j<nbj;j++ ) { 
            os << "("<<r[i][j].id<<","<<r[i][j].t<<")";
         }
         os<<std::endl;
      }
      os<<" delivery :"<<std::endl;
      for (unsigned int i=1;i<d.size();i++ ) {
         for (int j=0;j<nbj;j++ ) { 
            os << "("<<d[i][j].id<<","<<d[i][j].t<<")";
         }
         os<<std::endl;
      }
      os<< "Sum \n";
      for (unsigned int i=1;i<d.size();i++ ) {
            os << pk[i]<<std::endl;
      }
   }
   void Pack(Bob::Serialize &bs) const {
         Bob::Pack(bs,&r);
         Bob::Pack(bs,&d);
         Bob::Pack(bs,&pk);
         Bob::Pack(bs,&nbj);
   } 
   void UnPack(Bob::DeSerialize &ds) {
         Bob::UnPack(ds,&r);
         Bob::UnPack(ds,&d);
         Bob::UnPack(ds,&pk);
         Bob::UnPack(ds,&nbj);
   } 

};




