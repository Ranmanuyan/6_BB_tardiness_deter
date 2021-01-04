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
 *  File   : misc.h
 *  Author : B. Le Cun based on source of the ETSI Plugtest.
 *  Date   : Nov 11 2007.
 *  Comment: Misc Source file to solve the permutation flowshop using the Bob++ library
 */

#include <bobpp/bobpp>

#ifndef MISCHEADER
#define MISCHEADER

/** Template class which inherits std::vector
 * the goal is to add the pack and unpack methods
 */
template <class T> 
class pvector : public std::vector<T> {
public:
      pvector() : std::vector<T>() {}
      pvector(int s,T t=T()) : std::vector<T>(s,t) {}
      /// Pack method to serialize the BobNode
      virtual void Pack(Bob::Serialize &bs)  const {
         DBGAFF_ENV("pvector::Pack","Pack a pvector");
         const int s = std::vector<T>::size();
         Bob::Pack(bs,&s,1);
         for (int i=0;i<s;i++) {
            const T *v = &(*this)[i];
            Bob::Pack(bs,v);
         }
      }
      /// Unpack method to deserialize the BobNode
      virtual void UnPack(Bob::DeSerialize &bs)  {
         DBGAFF_ENV("pvector::UnPack","UnPack a pvector");
         int s;
         Bob::UnPack(bs,&s,1);
         std::vector<T>::resize(s);
         for (int i=0;i<s;i++) {
            T *v=&(*this)[i];
            Bob::UnPack(bs,v);
         }
      }
};

/** Permutation class used to store the branching information 
 */
class Permutation {
protected:
int nb;
int nbfree;
int r,d;
pvector<int> i2j;
pvector<int> j2i;
   public:
      /** Constructor
        */
      Permutation() : nb(0),nbfree(0),r(0),d(0),i2j(),j2i() {  }
      /** Constructor
        * @param s the size of the permutation
        */
      Permutation(int s) : nb(s),nbfree(s),r(0),d(s-1),i2j(s,-1),j2i(s,-1) { 
      }
      /** Constructor
        * @param p the permutation to copy
        */
      Permutation(const Permutation &p):nb(p.nb),nbfree(p.nbfree),r(p.r),d(p.d),i2j(p.i2j),j2i(p.j2i) { 
      }
      /** Destructor
        */
      ~Permutation() {  }
      /** Copy of of the permutation
        */
      virtual void copy(const Permutation &p) {
         nb = p.nb;
         nbfree = p.nbfree;
         r = p.r;
         d = p.d;
         i2j = p.i2j;
         j2i = p.j2i;
      }
      /** Method to initialize the permutation object
        * @param s the size of the permutation
        */
      void init(int s) {
         nbfree=nb=s;
         r=0;d=nb-1;
         i2j.resize(nb,-1); 
         j2i.resize(nb,-1); 
      }
      /** size return the size of the permutation
        */
      int size() const { return nb; }
      /** Fixe a variable say i next location 
        * @param i the value to put on the j location
        */
      void fixeR(int i) { fixe(i,r); r++; }
      /** Fixe a variable say i next location 
        * @param i the value to put on the j location
        */
      void fixeD(int i) { fixe(i,d); d--; }

      /** Fixe a variable say i to the place j
        * @param i the value to put on the j location
        * @param j the location
        */
      void fixe(int i,int j) { i2j[i]=j; j2i[j]=i; nbfree--;}
      /** get the value at the given location
        * @param i the location
        * @return the value of the location, -1 if the location is free
        */
      int geti(int j) const { return j2i[j]; }
      /** get the location a value
        * @param i the value 
        * @return the location of the value, -1 if the value is free
        */
      int getj(int i) const { return i2j[i]; }
      /** get the number of free locations/values
        * @return the number of free locations/values
        */
      int nbFree() const { return nbfree; }
      /** test if the value is fixed
        * @param i the value to test
        * @return true if the value is fixed
        */
      bool isfixed(int i) const { return i2j[i]!=-1; }
      /** test if the value is free
        * @param i the value to test
        * @return true if the value is free
        */
      bool isfree(int i) const { return i2j[i]==-1; }
      bool fixeOnR() const { return r<=(nb-1-d); }
      /** test there exists free locations/values
        * @return true if each value has been fixed
        */
      bool isallfixed() const { 
         return nbfree==0;
      }
      /// Display a Permutation
      virtual std::ostream &Prt(std::ostream &os) const {
        os << "values    :";
        for (int i=0;i<nb;i++ ) { os<<j2i[i]<<","; }
        os << std::endl;
        os << "locations :";
        for (int i=0;i<nb;i++ ) { os<<i2j[i]<<","; }
        os << std::endl;
        return os;
      }
      /// Pack method to serialize the BobNode
      virtual void Pack(Bob::Serialize &bs)  const {
         DBGAFF_ENV("Permutation::Pack","Pack a permutation");
         bs.Pack(&nb,1);
         bs.Pack(&nbfree,1);
         bs.Pack(&r,1);
         bs.Pack(&d,1);
         i2j.Pack(bs);
         j2i.Pack(bs);
      }
      /// Unpack method to deserialize the BobNode
      virtual void UnPack(Bob::DeSerialize &bs)  {
         DBGAFF_ENV("Permutation::UnPack","UnPack a permutation");
         bs.UnPack(&nb,1);
         bs.UnPack(&nbfree,1);
         bs.UnPack(&r,1);
         bs.UnPack(&d,1);
         i2j.UnPack(bs);
         j2i.UnPack(bs);
      }

};

#endif
