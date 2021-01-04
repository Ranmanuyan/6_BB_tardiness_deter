/*
 *        BOBO version 1.0:  Branch and Bound Optimization LiBrary
 *                    Opale Team of PRiSM laboratory
 *             University of Versailles St-Quentin en Yvelines.
 *                     Authors:  B. Le Cun F. Galea
 *                    (C) 2007 All Rights Reserved
 *
 *                              NOTICE
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted
 * provided that the above copyright notice appear in all copies and
 * that both the copyright notice and this permission notice appear in
 * supporting documentation.
 *
 * Neither the institutions (Versailles University, PRiSM Laboratory, 
 * the Opale Team), nor the Authors make any representations about the 
 * suitability of this software for any purpose.  This software is 
 * provided ``as is'' without express or implied warranty.
 *
 */

/* 
 *  File   : pfs.cpp
 *  Author : B. Le Cun based on source of the ETSI Plugtest.
 *  Date   : Nov 11 2007.
 *  Comment: Source file to solve the permutation flowshop using the Bob++ library
 */

#include"pfs.h"
#include"lb.h"



PFSInstance::PFSInstance() : Bob::BBInstance<PFSTrait>("pfs"),file(""),lobd(0),nbj(0),
         nbm(0),d(0), dueDate(0), lb1m(0){ 
      lobd=0;
}
   /// Constructor
PFSInstance::PFSInstance(std::string _file,int _lobd) : Bob::BBInstance<PFSTrait>("pfs"),
    file(_file),lobd(_lobd),nbj(0),nbm(0),d(0), dueDate(0){ 
   }
void PFSInstance::Init() { 
   read(); 
   lb1m= new OneMachine(this);  // must have this for mow. why?
}

// insert a root node and get the cost
PFSNode *PFSInstance::getSol() {
   PFSNode *sol=0;
   PFSNode *node;

   Bob::Permutation *per = new Bob::Permutation(); 
   per->init(nbj);
   
   for(int i=0; i<nbj; i++){
    per->fixe(i,i);
   }

   node = new PFSNode(this,*per);
   if ( sol==0 || sol->getEval()>node->getEval() ) {
    if ( sol!=0 ) delete sol; 
    sol=node;
    } else {
      delete node;
    }
   return sol;
}


void PFSInstance::read() {
   std::ifstream fr(file.c_str());

   if ( !fr.is_open() ) {
      std::cerr<< "Could not open file "<<file.c_str()<<" : exit\n";
      exit(1);
   }
   try {
     fr >> nbj;
     fr >> nbm;
     if ( nbj != 0 ) {
       d = new int *[ nbj ];
     }
     int ji;
     for ( int i = 0; i< nbj ; i++ ) {
       // fr >> ji;
       // fr >> ji;
       d[i]= new int[ nbm ];
       for ( int j= 0; j < nbm ; j++ ){
         fr >> ji;
         fr >> d[ i][ j];
       }
     }
     string due;
     fr >> due;
     if ( nbj != 0 ) {
       dueDate = new int [ nbj ];
     }
     for(int i=0; i< nbj; i++){
      fr >> ji;
      fr >> dueDate[i];
      fr >> ji;
      fr >> ji;
     }

   }
   catch( ... ) {
     cout << "Error while loading data from stream..." << endl;
   }
   for (int i=0;i<nbm;i++ ) {
     std::cout << "M["<<i<<"]:";
     for (int j=0;j<nbj;j++ ) {
       std::cout << d[j][i]<<",";
     }
     std::cout << std::endl;
   }
   std::cout << "Due Date: ";
   for (int j=0;j<nbj;j++ ) {
    std::cout << dueDate[j]<<",";
   }
   std::cout << std::endl;
   fr.close();

}

void PFSInstance::copy(const PFSInstance &pfi) {
   Bob::BBInstance<PFSTrait>::copy(pfi);
   file = pfi.file;
   lobd = pfi.lobd;
   nbj=pfi.nbj;nbm=pfi.nbm;
   d = new int *[ nbj ];
   dueDate = new int [nbj];
   for ( int i = 0; i< nbj ; i++ ) {
     d[i]= new int[ nbm ];
     dueDate[i] = pfi.dueDate[i];
     for ( int j= 0; j < nbm ; j++ )
       d[ i][ j] = pfi.d[i][j];
   }
   lb1m= new OneMachine(this);
}

   /// Pack method to serialize the BobNode
void PFSInstance::Pack(Bob::Serialize &bs) const  {
       DBGAFF_ENV("PFSInstance::Pack","Pack an instance");
	    Bob::BBInstance<PFSTrait>::Pack(bs);
       Bob::Pack(bs,&file);
       Bob::Pack(bs,&lobd);
       Bob::Pack(bs,&nbj);
       Bob::Pack(bs,&nbm);
       for (int i=0;i<nbj;i++ ) {
            Bob::Pack(bs,d[i],nbm);
            Bob::Pack(bs,&dueDate[i]);
       }
       Bob::Pack(bs,lb1m);
       DBGAFF_ENV("PFSInstance::Pack","------- Finish");
}
/// Unpack method to deserialize the BobNode
void PFSInstance::UnPack(Bob::DeSerialize &bs)  {
       DBGAFF_ENV("PFSInstance::UnPack","UnPack an instance");
	    Bob::BBInstance<PFSTrait>::UnPack(bs);
       Bob::UnPack(bs,&file);
       Bob::UnPack(bs,&lobd);
       Bob::UnPack(bs,&nbj);
       Bob::UnPack(bs,&nbm);
		 d = new int *[ nbj ];
     dueDate = new int [nbj];
       for (int i=0;i<nbj;i++ ) {
            d[i]=new int[ nbm ];
            Bob::UnPack(bs,d[i],nbm);
            Bob::UnPack(bs,&dueDate[i]);
       }
       lb1m=new  OneMachine();
       Bob::UnPack(bs,lb1m);
       DBGAFF_ENV("PFSInstance::UnPack","------- Finish");
}

void Scheduled::fus() {
   int max=0;
   //Prt(std::cout);
   for (unsigned int i=0;i<r.size();i++ ) {
      putr(i,getr(i)+getd(i));
      //std::cout<<"m:"<<i<<":"<<getr(i)<<"max:"<<max<<std::endl;
      if ( getr(i)>max ) {
         max=getr(i);
      }
      putd(i,max);
      putr(i,max);
   }
}

int Scheduled::getCost(const PFSInstance &pfi, const Bob::Permutation &per){
  
  vector<int> j2i =  per.get_j2i();
  vector<int> PS = vector<int>(j2i.begin(), j2i.begin()+per.size()-per.nbFree());
  const PFSInstance *pi = &pfi;
  int totalTardiness = calTardiness(PS, pi);


  j2i.clear();
  PS.clear();
  
  // delete pi;

  return totalTardiness;
}

void Scheduled::FixeR(const PFSInstance &pfi, int j) {
   putr(0,getr(0)+pfi.get(0,j));
   for (unsigned int i=1;i<r.size();i++ ) {
      putr(i,max(getr(i),getr(i-1))+pfi.get(i,j));
   }
}
// void Scheduled::FixeD(const PFSInstance &pfi, int j) {
//    int s = d.size()-1;
//    putd(s,getd(s)+pfi.get(s,j));
//    for (int i=s-1;i>-1;i-- ) {
//       putd(i,max(getd(i),getd(i+1))+pfi.get(i,j));
//    }
// }

PFSNode::PFSNode(const PFSInstance *pi,const Bob::Permutation &per):Bob::BBIntMinNode(0),nbj(pi->nbj),nbm(pi->nbm),per(nbj),sc(nbm) { 
   for (int j=0;j<nbj;j++ ) {
      FixeR(per.geti(j),*pi);
      dist()++;
   }
   const PFSInstance &pfi = *pi;
   setEval(sc.getCost(pfi, per));

}

void PFSNode::FixeR(int j,const PFSInstance &pfi) {
   per.fixeR(j);
   sc.FixeR(pfi,j);
   if ( per.isallfixed() ) {
      sc.fus();
      setEval(sc.getCost(pfi, per));
      setSol();
   }
}

// void PFSNode::FixeD(int j,const PFSInstance &pfi) {
//    per.fixeD(j);
//    sc.FixeD(pfi,j);
//    if ( per.isallfixed() ) {
//       sc.fus();
//       setEval(sc.getCost());
//       setSol();
//    }
// }

void PFSNode::eval(const PFSInstance *pfi,int minCost) {
   // int lb5;
   // setEval(lb1);
   if(minCost ==-1){
    setEval(0);
   }else{
    int lb1 = pfi->lb1m->LB2(pfi,per,sc);
    setEval(lb1);
   }

   
   // if ( minCost==-1 || minCost>lb1 ) {
   //    lb5 = pfi->lb2m->LB5(pfi,per,sc);
   //    setEval(lb5);
   // } 
}

#if 1
bool PFSGenChild::operator()(PFSNode *p) {
   std::cout << "------Debut GenChild---------------------------------------------\n";
   p->Prt(std::cout); // pfs.h row:189
   for (int j=0;j<inst->nbj;j++) {
      if ( p->perm().isfree(j) )  {         
         std::cout << "------Fils "<<j<<"--------------------------------------\n";
         PFSNode *nf=new PFSNode(p);
         nf->FixeR(j,*inst);
         if ( nf->isSol() ) {
            std::cout << "***********Solution :"<<nf->getEval()<<std::endl;
            nf->Prt(std::cout);
         } else {
            nf->eval(inst,algo->getGoal()->getBest());
            std::cout << " Generate:"<<j<<" on "<<nf->dist()<<std::endl; // nf->dist is location index
            nf->Prt(std::cout);
            nf->dist()++;
         }
         algo->Search(nf);
      } 
   }
   std::cout << "------Fin GenChild---------------------------------------------\n";
}
#else
bool PFSGenChild::operator()(PFSNode *p) {
   bool fl=false;
   PFSNode *nf;
   AFFMESS(std::cout << "------Debut GenChild---------------------------------------------\n");
   AFFMESS(p->Prt(std::cout));
   if ( p->perm().nbFree()==2 ) {
      int i1=-1,i2=-1;
      for (int j=0;j<inst->nbj;j++ ) {
         if ( p->perm().isfree(j) ) {
            if ( i1==-1 )  i1=j;
            else if ( i1!=-1 && i2==-1 )  i2=j;
            else std::cout <<"strange\n";
         }
      }
      nf=new PFSNode(p);
      nf->FixeR(i1,*inst);
      nf->FixeR(i2,*inst);
      AFFMESS(std::cout << "***********Solution :"<<nf->getEval()<<std::endl);
      //nf->Prt(std::cout);
      fl += Search(nf);
      nf=new PFSNode(p);
      nf->FixeR(i2,*inst);
      nf->FixeR(i1,*inst);
      AFFMESS(std::cout << "***********Solution :"<<nf->getEval()<<std::endl);
      //nf->Prt(std::cout);
      fl +=Search(nf) || fl;
      return fl;
   }
   if ( p->perm().fixeOnR() ) {
      for (int j=0;j<inst->nbj;j++ ) {
         if ( p->perm().isfree(j) )  {
            AFFMESS(std::cout << "------Fils R"<<j<<"--------------------------------------\n");
            int ir = p->dist()+1-Bob::core::get_depth_rec();
            if ( ir>=0 ) {
              nf = get_rec_node(ir);
              if ( nf==0 ) {
                 nf = new PFSNode(p); 
                 declare_as_rec(nf,ir);
              }  
              nf->copy(p);
              nf->addRef();
            } else {
              nf = new PFSNode(p); 
            }
            nf->FixeR(j,*inst);
            algo->start_eval(nf,p);

            nf->eval(inst,algo->getGoal()->getBest());
            algo->end_eval(nf);
            AFFMESS(std::cout << " Generate:"<<j<<" on IN"<<std::endl);
            AFFMESS(nf->Prt(std::cout));
            nf->dist()++;
            fl +=Search(nf);
         }
      } 
   } else {
      for (int j=0;j<inst->nbj;j++ ) {
         if ( p->perm().isfree(j) )  {
            AFFMESS(std::cout << "------Fils D"<<j<<"--------------------------------------\n");
            nf=new PFSNode(p);
            nf->FixeD(j,*inst);
            algo->start_eval(nf,p);
            nf->eval(inst,algo->getGoal()->getBest());
            algo->end_eval(nf);
            AFFMESS(std::cout << " Generate:"<<j<<" on OUT"<<std::endl);
            AFFMESS(nf->Prt(std::cout));
            nf->dist()++;
            fl +=Search(nf);
         } 
      }
   }
   AFFMESS(std::cout << "------Fin GenChild---------------------------------------------\n");
   return fl;
}

#endif




int main(int n, char ** v) {
  Bob::core::opt().add(std::string("--pfs"), Bob::Property("-f", "Problem File", "../taillard/DD_Ta000.txt"));
  Bob::core::opt().add(std::string("--pfs"), Bob::Property("-lb", "Lower bound", 0));
#ifdef Atha
  Bob::AthaBBAlgoEnvProg<PFSTrait> env;
  Bob::AthaEnvProg::init(n, v);
  Bob::core::Config(n, v);
#elif defined(Threaded)
  Bob::ThrBBAlgoEnvProg<PFSTrait> env;
  Bob::ThrEnvProg::init();
  Bob::core::Config(n, v);
  Bob::ThrEnvProg::start();
#elif defined(MPxMPI)
  Bob::mpx::EnvProg::init(n, v);
  Bob::mpx::BBAlgoEnvProg<PFSTrait,Bob::mpx::lb_strat_2t,Bob::mpx::topo_a2a,long> env;
  Bob::core::Config(n, v);
  Bob::mpx::EnvProg::start();
#else
  Bob::SeqBBAlgoEnvProg<PFSTrait> env;
  Bob::core::Config(n, v);
#endif

  PFSInstance *Instance = new PFSInstance(Bob::core::opt().SVal("--pfs", "-f"),Bob::core::opt().NVal("--pfs", "-lb"));
  // Initialize the instance.

  // Solve the instance i.e. all the tree
  std::cout << "-------- Start to solve the entire tree" << std::endl;
  env(Instance);
 
#ifdef Atha
  Bob::AthaEnvProg::end();
#elif defined(Threaded)
  Bob::ThrEnvProg::end();
#elif defined(MPxMPI)
  Bob::mpx::EnvProg::end();
#endif
  delete Instance;
  Bob::core::End();
  exit(0);
}

