#include "Cryptosystem.h"

#define N 8
#define Q 97/*195161//7919//1049//17*/

using namespace std;
using namespace NTL;

void Add(ZZX& sum, ZZX& a, ZZX& b, ZZ q) {
  sum.SetLength(N);
  for(int i = 0; i < N; i++) {
    sum[i] = (a[i] + b[i]) % q;
    if(sum[i] > ((q-1)/2)) sum[i] = sum[i] - q;
  }
}

void MyProd(ZZX& prod, ZZX& a, ZZX& b, ZZ q) {
  prod.SetLength(N);
  for(int i = 0; i < N; i++)
    prod[i] = (a[i] & b[i]) % q;
}

void Sub(ZZX& sub, ZZX& a, ZZX& b, ZZ q) {
  sub.SetLength(N);
  for(int i = 0; i < N; i++)
    sub[i] = (a[i] - b[i]) % q;
}

void ProdByScalar(ZZX& prod, ZZX& a, ZZ b, ZZX& f, ZZ Zq) {
  prod.SetLength(N);

  prod = a * b;

  prod = prod % f;

  for(int i = 0; i < N; i++){
    prod[i] = (prod[i] % Zq) % Zq;
    if(prod[i] > ((Zq-1)/2)) prod[i] = prod[i] - Zq;
  }
}

void Prod(ZZX& prod, ZZX& a, ZZX& b, ZZX& f, ZZ Zq) {
  prod.SetLength(N);
  NTL::MulMod(prod, a, b, f);

  for(int i = 0; i < N; i++){
    prod[i] = prod[i] % Zq;
    if(prod[i] > ((Zq-1)/2)) prod[i] = prod[i] - Zq;
  }
}

void Encode(ZZX& m_bar, ZZX& m, ZZ Zq) {
  m_bar.SetLength(N);

  for(int i = 0; i < N; i++)
    m_bar[i] = m[i] > 0 ? (to_ZZ((Zq - 1) / 2)) : to_ZZ(0);
}

void Decode(ZZX& m_d, ZZX& m_prime, ZZ Zq) {
  m_d.SetLength(N);
  ZZ lbound = to_ZZ((Zq-1)/4);
  ZZ ubound = to_ZZ(3*lbound);

  for(int i = 0; i < N; i++) {
    cout << "m_prime[i]: " << m_prime[i] << ", lbound: " << lbound << ", ubound: " << ubound << endl;
      if(m_prime[i] >= lbound && m_prime[i] < ubound)
          m_d[i] = (int32_t) 1;
      else
          m_d[i] = (int32_t) 0;
  }
}

int main(int argc, char** argv) {
  //ZZ_p::init(to_ZZ(97));
  ZZ q = to_ZZ(97);
  ZZX s, e, a1, a0, f, m;

  s.SetLength(8);
  e.SetLength(8);
  a1.SetLength(8);
  a0.SetLength(8);
  f.SetLength(9);
  m.SetLength(8);

  f[8] = 1; f[0] = 1;
  s[7] = -1;s[6] = -1;s[5] = -1;s[4] = 1;s[3] = 1;s[2] = 1;s[1] = 1;s[0] = -1;
  e[7] = 0;e[6] = -2;e[5] = 0;e[4] = 0;e[3] = -2;e[2] = 2;e[1] = -1;e[0] = 1;
  a1[7] = 27;a1[6] = -11;a1[5] = -33;a1[4] = 41;a1[3] = -18;a1[2] = -5;a1[1] = -37;a1[0] = -16;
  m[7] = 0;m[6] = 1;m[5] = 0;m[4] = 1;m[3] = 1;m[2] = 0;m[1] = 0;m[0] = 1;

  ZZX a1s, e2;
  a1s.SetLength(8);
  e2.SetLength(8);

  Prod(a1s, a1, s, f, q);
  ProdByScalar(e2, e, to_ZZ(2), f, q);
  Add(a0, a1s, e2, q);
  ProdByScalar(a0, a0, to_ZZ(-1), f, q);

  cout << "a0: " << a0 << endl;

  ZZX c0, c1, u, g, h;

  c0.SetLength(8);
  c1.SetLength(8);
  u.SetLength(8);
  g.SetLength(8);
  h.SetLength(8);

  u[7] = 0;u[6] = -2;u[5] = 3;u[4] = 0;u[3] = 2;u[2] = 0;u[1] = -1;u[0] = 0;
  g[7] = 0;g[6] = -1;g[5] = 0;g[4] = 0;g[3] = 0;g[2] = -1;g[1] = 2;g[0] = 0;
  h[7] = -1;h[6] = 0;h[5] = 1;h[4] = 1;h[3] = 0;h[2] = 0;h[1] = 1;h[0] = 1;

  ZZX a0u, g2;

  //Prod(a0u, a0, u, f, q);
  //ProdByScalar(g2, g, to_ZZ(2), f, q);
  //Add(c0, a0u, g2, q);
  //Add(c0, c0, m, q);
  c0 = ((a0 * u) % f) + ((2 * g) % f) + m;
  for(int i = 0; i < N; i++){
    c0[i] = c0[i] % q;
    if(c0[i] > ((q-1)/2)) c0[i] = c0[i] - q;
  }

  cout << "c0: " << c0 << endl;

  /*c1 = ((a1 * u) % f) + ((2 * h) % f);

  cout << "c1: " << c1 << endl;

  ZZ_pX m_tilde;
  m_tilde.SetLength(8);

  m_tilde = c0 + ((c1 * s) % f);

  cout << "m_tilde: " << m_tilde << endl;*/

  return 0;
}

int main1(int argc, char** argv) {
  HElium::Cryptosystem cs;
  cs.KeyGen(N, Q);
  ZZX m1, m2, md, mm, ct1[2], ct2[2], ctAdd[2], ctProd[3];

  m1.SetLength(N);
  //m1[3] = to_ZZ(1); m1[2] = to_ZZ(0); m1[1] = to_ZZ(0); m1[0] = to_ZZ(1);
  m1[3] = 0; m1[2] = 0; m1[1] = 0; m1[0] = 1;
  cout << "m1: " << m1 << endl;

  m2.SetLength(N);
  //m2[3] = to_ZZ(1); m2[2] = to_ZZ(0); m2[1] = to_ZZ(1); m2[0] = to_ZZ(0);
  m2[3] = 0; m2[2] = 0; m2[1] = 0; m2[0] = 1;
  cout << "m2: " << m2 << endl;

  cs.Encrypt(ct1, m1);
  //cs.Encrypt(ct2, m2);

  //cs.HEAdd(ctAdd, ct1, ct2);

  //cs.Decrypt(md, ctAdd);

  //cs.NOT(ct1);

  cs.Decrypt(mm, ct1);

  //cs.AND(ctAdd, ct1, ct2);

  //cs.Decrypt(mm, ctAdd);

  //cs.HEMult(ctProd, ct1, ct2);
  //cs.Decrypt3(mm, ctProd);

  cout << "= END =================================================================" << endl;
  return 0;
}

int main0(int argc, char** argv) {
  ZZX f, r1, r2, a, a_r2, p, m1, m2, m1_bar, e1, e2, e3, cm11, cm12, a_e1, p_e1, m_prime, cm11_r2, m_d;
  ZZX m2_bar, cm21, cm22, m2_prime, cm21_r2, m2_d, sumc1, sumc2;
  f.SetLength(N + 1);
  f[4] = to_ZZ(1);
  f[0] = to_ZZ(1);

  cout << "n: " << N << endl;
  cout << "q: " << Q << endl;
  cout << "f(x): " << f << endl;



  cout << "= KeyGen: =============================================================" << endl;

  r1.SetLength(N);
  r1[3] = to_ZZ(11); r1[2] = to_ZZ(3); r1[1] = to_ZZ(2); r1[0] = to_ZZ(6);
  cout << "r1: " << r1 << endl;

  r2.SetLength(N);
  r2[3] = to_ZZ(2); r2[2] = to_ZZ(5); r2[1] = to_ZZ(8); r2[0] = to_ZZ(7);
  cout << "r2: " << r2 << endl;

  a.SetLength(N);
  //a[3] = to_ZZ(3); a[2] = to_ZZ(1); a[1] = to_ZZ(11); a[0] = to_ZZ(8);
  cout << "a: " << a << endl;

  Prod(a_r2, a, r2, f, to_ZZ(Q));

  cout << "a_r2: " << a_r2 << endl;

  Sub(p, r1, a_r2, to_ZZ(Q));

  cout << "p: " << p << endl;

  cout << "= Encrypt: ============================================================" << endl;

  m1.SetLength(N);
  //m1[3] = to_ZZ(1); m1[2] = to_ZZ(0); m1[1] = to_ZZ(0); m1[0] = to_ZZ(1);
  //m1[3] = 0; m1[2] = 1; m1[1] = 1; m1[0] = 1;
  m1[3] = 0; m1[2] = 0; m1[1] = 0; m1[0] = 1;
  cout << "m1: " << m1 << endl;

  m2.SetLength(N);
  //m2[3] = to_ZZ(1); m2[2] = to_ZZ(0); m2[1] = to_ZZ(1); m2[0] = to_ZZ(0);
  //m2[3] = 0; m2[2] = 0; m2[1] = 1; m2[0] = 1;
  m2[3] = 0; m2[2] = 0; m2[1] = 0; m2[0] = 1;
  cout << "m2: " << m2 << endl;

  Encode(m1_bar, m1, to_ZZ(Q));
  cout << "m1_bar: " << m1_bar << endl;

  e1.SetLength(N);
  e1[3] = to_ZZ(1); e1[2] = to_ZZ(3); e1[1] = to_ZZ(1); e1[0] = to_ZZ(6);
  //e1[3] = 0; e1[2] = 0; e1[1] = 0; e1[0] = 0;
  cout << "e1: " << e1 << endl;

  e2.SetLength(N);
  e2[3] = to_ZZ(4); e2[2] = to_ZZ(2); e2[1] = to_ZZ(3); e2[0] = to_ZZ(1);
  //e2[3] = 0; e2[2] = 0; e2[1] = 0; e2[0] = 0;
  cout << "e2: " << e2 << endl;

  e3.SetLength(N);
  e3[3] = to_ZZ(1); e3[2] = to_ZZ(4); e3[1] = to_ZZ(2); e3[0] = to_ZZ(1);
  //e3[3] = 0; e3[2] = 0; e3[1] = 0; e3[0] = 0;

  cout << "e3: " << e3 << endl;

  Prod(a_e1, a, e1, f, to_ZZ(Q));
  cout << "a_e1: " << a_e1 << endl;

  Add(cm11, a_e1, e2, to_ZZ(Q));
  cout << "cm11: " << cm11 << endl;

  Prod(p_e1, p, e1, f, to_ZZ(Q));
  cout << "p_e1: " << p_e1 << endl;

  Add(cm12, p_e1, e3, to_ZZ(Q));
  cout << "cm12: " << cm12 << endl;
  Add(cm12, cm12, m1_bar, to_ZZ(Q));
  cout << "cm12: " << cm12 << endl;

/////////////////////////////////////////////////////////////////////////////////

  Encode(m2_bar, m2, to_ZZ(Q));
  cout << "m2_bar: " << m2_bar << endl;

  Add(cm21, a_e1, e2, to_ZZ(Q));
  cout << "cm21: " << cm21 << endl;

  Add(cm22, p_e1, e3, to_ZZ(Q));
  cout << "cm22: " << cm22 << endl;
  Add(cm22, cm22, m2_bar, to_ZZ(Q));
  cout << "cm22: " << cm22 << endl;

  cout << "= Decrypt: ============================================================" << endl;

  Prod(cm11_r2, cm11, r2, f, to_ZZ(Q));
  cout << "cm11_r2: " << cm11_r2 << endl;

  Add(m_prime, cm11_r2, cm12, to_ZZ(Q));
  cout << "m_prime: " << m_prime << endl;

  Decode(m_d, m_prime, to_ZZ(Q));
  cout << "m_d: " << m_d << endl;

/////////////////////////////////////////////////////////////////////////////////

  Prod(cm21_r2, cm21, r2, f, to_ZZ(Q));
  cout << "cm21_r2: " << cm21_r2 << endl;

  Add(m2_prime, cm21_r2, cm22, to_ZZ(Q));
  cout << "m2_prime: " << m2_prime << endl;

  Decode(m2_d, m2_prime, to_ZZ(Q));
  cout << "m2_d: " << m2_d << endl;

  cout << "= Sum: ================================================================" << endl;

  Add(sumc1, cm11, cm21, to_ZZ(Q));
  Add(sumc2, cm12, cm22, to_ZZ(Q));

  ZZX sumc1_r2, ms_prime, ms_d;

  Prod(sumc1_r2, sumc1, r2, f, to_ZZ(Q));
  cout << "sumc1_r2: " << sumc1_r2 << endl;

  Add(ms_prime, sumc1_r2, sumc2, to_ZZ(Q));
  cout << "ms_prime: " << ms_prime << endl;

  Decode(ms_d, ms_prime, to_ZZ(Q));
  cout << "ms_d: " << ms_d << endl;

  cout << "= Mult: ===============================================================" << endl;

  ZZX mulc1, mulc2, mulc3, mulc2p1, mulc2p2;

  Prod(mulc1, cm11, cm21, f, to_ZZ(Q));
  cout << "mulc1: " << mulc1 << endl;

  Prod(mulc2p1, cm11, cm22, f, to_ZZ(Q));
  cout << "mulc2p1: " << mulc2p1 << endl;

  Prod(mulc2p2, cm21, cm12, f, to_ZZ(Q));
  cout << "mulc2p2: " << mulc2p2 << endl;

  Add(mulc2, mulc2p1, mulc2p2, to_ZZ(Q));
  cout << "mulc2: " << mulc2 << endl;

  Prod(mulc3, cm12, cm22, f, to_ZZ(Q));
  cout << "mulc3: " << mulc3 << endl;

  ZZX r2_2, mulc1_r2_2, mulc2_r2, mm_prime, mm_d;

  Prod(r2_2, r2, r2, f, to_ZZ(Q));
  cout << "r2_2: " << r2_2 << endl;

  Prod(mulc1_r2_2, mulc1, r2_2, f, to_ZZ(Q));
  cout << "mulc1_r2_2: " << mulc1_r2_2 << endl;

  Prod(mulc2_r2, mulc2, r2, f, to_ZZ(Q));
  cout << "mulc2_r2: " << mulc2_r2 << endl;

  Add(mm_prime, mulc1_r2_2, mulc2_r2, to_ZZ(Q));
  cout << "mm_prime: " << mm_prime << endl;
  Add(mm_prime, mm_prime, mulc3, to_ZZ(Q));
  cout << "mm_prime: " << mm_prime << endl;

  Decode(mm_d, mm_prime, to_ZZ(Q));
  cout << "mm_d: " << mm_d << endl;

cout << "= MY SH!T: ===============================================================" << endl;

MyProd(sumc1, cm11, cm21, to_ZZ(Q));
MyProd(sumc2, cm12, cm22, to_ZZ(Q));

//ZZX sumc1_r2, ms_prime, ms_d;

Prod(sumc1_r2, sumc1, r2, f, to_ZZ(Q));
cout << "myprodc1_r2: " << sumc1_r2 << endl;

Add(ms_prime, sumc1_r2, sumc2, to_ZZ(Q));
cout << "myprodms_prime: " << ms_prime << endl;

Decode(ms_d, ms_prime, to_ZZ(Q));
cout << "myprod: " << ms_d << endl;

  return 0;
}
