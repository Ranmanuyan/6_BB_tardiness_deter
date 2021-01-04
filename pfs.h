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
 *  File   : flowshop.cpp
 *  Author : B. Le Cun based on source of the ETSI Plugtest.
 *  Date   : Nov 11 2007.
 *  Comment: Source file to solve the permutation flowshop using the Bob++ library
 */


#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>

#include <bobpp/bobpp>
#ifdef Threaded
#  include <bobpp/thr/thr.h>
#elif defined(Atha)
#  include <bobpp/atha/atha.h>
#elif defined(MPxMPI)
#  include <bobpp/mpx/mpx.hpp>
#endif



#ifndef PFSHEADER
#define PFSHEADER

/** total tardiness minimization of flow shop scheduling
    @author Lei Liu
    @version 1.0
    17 June 2020
**/


class PFSInstance;
class PFSNode;
class PFSGenChild;
class PFSStat;

class PFSTrait {
  public:    
    typedef PFSNode Node;
    typedef PFSInstance Instance;
    typedef PFSGenChild GenChild;
    typedef Bob::BBAlgo<PFSTrait> Algo;
    typedef Bob::BBGoalBest<PFSTrait> Goal;
    typedef Bob::DepthEPri<PFSNode> PriComp;
    typedef Bob::BBStat Stat;
};

#define AFFMESS(m) 
//#define AFFMESS(m) m

/** Class used to compute the cost of a partial solution.
 */
class Scheduled {
Bob::pvector<int> r;
Bob::pvector<int> d;

public:
      Scheduled() : r(),d() {} 
      Scheduled(int n) : r(n,0),d(n,0) { } 
      Scheduled(const Scheduled &md) : r(md.r),d(md.d) {} 
      virtual ~Scheduled() {} 
      virtual void copy(const Scheduled &md){
         r = md.r; d=md.d;
      }
      int getr(int i) const { return r[i]; }
      void putr(int i,int v) { r[i]=v; }
      int getd(int i) const { return d[i]; }
      void putd(int i,int v) { d[i]=v; }
      void fus();
      
      //the cost of a partial schedule
      int getCost(const PFSInstance &pfi, const Bob::Permutation &per); //return getr(r.size()-1);
      

      /// Pack method to serialize the BobNode
      virtual void Pack(Bob::Serialize &bs)  const {
         DBGAFF_ENV("Scheduled::Pack","Pack a Scheduled");
         r.Pack(bs);
         d.Pack(bs);
      }
      /// Unpack method to deserialize the BobNode
      virtual void UnPack(Bob::DeSerialize &bs)  {
         DBGAFF_ENV("Scheduled::UnPack","UnPack a Scheduled");
         r.UnPack(bs);
         d.UnPack(bs);
      }
      void FixeR(const PFSInstance &pfi, int i);
      // void FixeD(const PFSInstance &pfi, int i);
      /// Display the Machine Date
      // virtual ostream &Prt(ostream &os) const {
      //   unsigned int i;
      //   os << "Cr:"; for (i=0;i<r.size();i++ ) { os<<getr(i)<<","; } os << std::endl;
      //   os << "Cd:"; for (i=0;i<r.size();i++ ) { os<<getd(i)<<","; } os << std::endl;
      //   return os;
      // }

};

/** the Node of the search
 */
class PFSNode : public Bob::BBIntMinNode {
protected:
int nbj;
int nbm;
Bob::Permutation per;
Scheduled sc;
public:
    /// Contructor
    PFSNode() : Bob::BBIntMinNode(),nbj(0),nbm(0),per(),sc() { }
    /// Contructor
    PFSNode(int _nbj,int _nbm) : Bob::BBIntMinNode(0),nbj(_nbj),nbm(_nbm),per(_nbj),sc(_nbm) { 
    }
    /// Contructor
    PFSNode(int _nbj,int _nbm,int cost) : Bob::BBIntMinNode(cost),nbj(_nbj),nbm(_nbm),per(_nbj),sc(_nbm) { 
    }
    PFSNode(const PFSNode &p) : Bob::BBIntMinNode(p),nbj(p.nbj),nbm(p.nbm),per(p.per),sc(p.sc) {
    }
    /// Contructor
    PFSNode(PFSNode *p) : Bob::BBIntMinNode(*p),nbj(p->nbj),nbm(p->nbm),per(p->per),sc(p->sc) { 
    }
    /// Contructor
    PFSNode(const PFSInstance *pi,const Bob::Permutation &per);
    /// Fixe the job i on the next free location
    void FixeR(int i,const PFSInstance &pfi);
    /// Fixe the job i on the last free location
    // void FixeD(int i,const PFSInstance &pfi);

    /// Destructor
    virtual ~PFSNode() {
    }
    virtual void copy(const PFSNode &pn) {
      Bob::BBIntMinNode::copy(pn);
      nbj= pn.nbj;
      nbm= pn.nbm;
      per.copy(pn.per);
      sc.copy(pn.sc);
    }
    /// get the permutation
    Bob::Permutation &perm() { return per; }
    /// get the Scheduled 
    Scheduled &sched() { return sc; }
    /// Pack method to serialize the BobNode
    virtual void Pack(Bob::Serialize &bs)  const {
       DBGAFF_ENV("PFSNode::Pack","Pack a node");
	    Bob::BBIntMinNode::Pack(bs);
       Bob::Pack(bs,&nbj,1);
       Bob::Pack(bs,&nbm,1);
       Bob::Pack(bs,&per);
       Bob::Pack(bs,&sc);
    }
    /// Unpack method to deserialize the BobNode
    virtual void UnPack(Bob::DeSerialize &bs)  {
       DBGAFF_ENV("PFSNode::UnPack","UnPack a node");
	    Bob::BBIntMinNode::UnPack(bs);
       Bob::UnPack(bs,&nbj,1);
       Bob::UnPack(bs,&nbm,1);
       Bob::UnPack(bs,&per);
       Bob::UnPack(bs,&sc);
    }

    /// evaluate the node
    void eval(const PFSInstance *pfi,int minCost);
    /// Display a PFS Node
    // display the locations and the cost for partial schedule and full schedule
    virtual ostream &Prt(ostream &os) const {
        if ( isPb() ) {
           os << "["<<getEval()<<"]\n ";
           per.Prt(os);
           // sc.Prt(os);
        } else {
           os << "Solution Cost :"<<getEval()<<"\n";
           per.Prt(os);
        }
        return os;
    }
    /** Fill the strbuff with node information
      * @param st the strbuff to fill.
      */
    virtual void log_space(Bob::strbuff<> &st) {
      Bob::BBIntMinNode::log_space(st);
      if ( dist()==0 ) { st<<",\"root\""; }
      else {
        if ( perm().fixeOnR() ) 
           st <<",\""<<per.getFacD()<<","<<per.getLocD()<<"\"";
        else 
           st <<",\""<<per.getFacR()<<","<<per.getLocR()<<"\"";
           
      }
    }

};

/** the genchild class for the permutation flowshop
 */
class PFSGenChild : public Bob::BBGenChild<PFSTrait>  {
        public:
          PFSGenChild(const PFSInstance *_inst, PFSTrait::Algo *_sch) : 
                 Bob::BBGenChild<PFSTrait>(_inst,_sch) {}
          virtual ~PFSGenChild() {}
          virtual bool operator()(PFSNode *n);
};



class OneMachine;

/// The Instance class for the PFS.
class PFSInstance : public Bob::BBInstance<PFSTrait> {
public:
/// The file name of the problem instance
std::string file;
int lobd;
int nbj,nbm;
int ** d;
int * dueDate;
OneMachine *lb1m;  //  instance Lower bound


   /// Constructor
   PFSInstance() ;
   /// Constructor
   PFSInstance(std::string _file,int _lobd);
   /// Destructor
   virtual ~PFSInstance() {}
   virtual void copy(const PFSInstance &pfi);
   /// Initialization of the instance
   void Init();
   /// the method that initialize and insert the root nodes.
   void InitAlgo(PFSTrait::Algo *al) const {
       PFSNode *qn;
       qn = new PFSNode(nbj,nbm);
       qn->eval(this,-1);  //this is the LB limit, give it as small as possible based on your objective function
       std::cout << "Root Evaluation :"<<qn->getEval()<<std::endl;
       //al->log_space(qn);
       al->GPQIns(qn);
   }
   /// the method that initialize and insert the initial upper bound.
   PFSNode *getSol();
   /// reads the data from files.
   void read();
   unsigned long get(int i,int j) const { return d[j][i]; }
   /// Pack method to serialize the BobNode
   virtual void Pack(Bob::Serialize &bs) const;
   /// Unpack method to deserialize the BobNode
   virtual void UnPack(Bob::DeSerialize &bs);
};

#define max(a,b) ((a)<(b)?(b):(a))

#endif
