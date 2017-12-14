#include <iostream>
#include <time.h>
#include <random>

#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>

#ifndef CRYPTOSYSTEM_H
#define CRYPTOSYSTEM_H

using namespace std;
using namespace NTL;

namespace HElium {

class Cryptosystem {
public:
  ZZ_pX pk[2];
  int32_t n;
  long q;
  //Cryptosystem(int32_t n, long q, bool generateKeys);
  Cryptosystem();
  void KeyGen(int32_t n, long q);
  void Encrypt(ZZ_pX ct[],ZZ_pX& m);
  void Decrypt(ZZ_pX&m, ZZ_pX ct[]);
  void Decrypt3(ZZ_pX&m, ZZ_pX ct[]);
  void HEAdd(ZZ_pX sum[], ZZ_pX ct1[], ZZ_pX ct2[]);
  void HEMult(ZZ_pX prod[], ZZ_pX ct1[], ZZ_pX ct2[]);
  void AND(ZZ_pX sum[], ZZ_pX ct1[], ZZ_pX ct2[]);
  void XOR();
  void NAND();
  void NOT(ZZ_pX ct[]);
protected:
  ZZ_pX f;
  void PolyAdd(ZZX& sum, ZZX& a, ZZX& b);
  void PolySub(ZZX& sub, ZZX& a, ZZX& b);
  void PolyProdByScalar(ZZX& prod, ZZX& a, ZZ b);
  void PolyProd(ZZX& prod, ZZX& a, ZZX& b);
  void PolyRandom(ZZX& poly, long mod);
  void PolyRandom(ZZX& poly);
private:
  ZZ_pX sk;
  void Encode(ZZX& m_bar, ZZX& m);
  void Decode(ZZX& m_d, ZZX& m_prime);
};

}

#endif
