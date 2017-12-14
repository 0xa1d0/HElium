#include "Cryptosystem.h"

using namespace std;
using namespace NTL;
using namespace HElium;

void Cryptosystem::KeyGen(int32_t n, long q) {
  this->n = n;
  this->q = q;

  ZZ_p::init(to_ZZ(this->q));

  this->f.SetLength(this->n + 1);
  this->f[this->n] = 1;
  this->f[0] = 1;

  ZZ_pX s, e, a1, a0;

  cout << "= KeyGen: =============================================================" << endl;

  a1.SetLength(this->n);
  a0.SetLength(this->n);
  e.SetLength(this->n);
  s.SetLength(this->n);

  s[7] = -1;s[6] = -1;s[5] = -1;s[4] = 1;s[3] = 1;s[2] = 1;s[1] = 1;s[0] = -1;
  e[7] = 0;e[6] = -2;e[5] = 0;e[4] = 0;e[3] = -2;e[2] = 2;e[1] = -1;e[0] = 1;
  a1[7] = 27;a1[6] = -11;a1[5] = -33;a1[4] = 41;a1[3] = -18;a1[2] = -5;a1[1] = -37;a1[0] = -16;

  a0 = -1 * (((a1 * s) % this->f) + ((2 * e) % this->f));

  this->sk = s;
  cout << "sk: " << s << endl;
  this->pk[0] = a0;
  this->pk[1] = a1;
  cout << "pk: {" << a0 << ", " << a1 << "}" << endl;
}

void Cryptosystem::Encrypt(ZZ_pX ct[],ZZ_pX& m) {
  ZZ_pX c0, c1, u, g, h;

  c0.SetLength(this->n);
  c1.SetLength(this->n);
  u.SetLength(this->n);
  g.SetLength(this->n);
  h.SetLength(this->n);

  u[7] = 0;u[6] = -2;u[5] = 3;u[4] = 0;u[3] = 2;u[2] = 0;u[1] = -1;u[0] = 0;
  g[7] = 0;g[6] = -1;g[5] = 0;g[4] = 0;g[3] = 0;g[2] = -1;g[1] = 2;g[0] = 0;
  h[7] = -1;h[6] = 0;h[5] = 1;h[4] = 1;h[3] = 0;h[2] = 0;h[1] = 1;h[0] = 1;

  c0 = ((this->pk[0] * u) % f) + ((2 * g) % this->f) + m;

  cout << "c0: " << c0 << endl;

  c1 = ((this->pk[1] * u) % f) + ((2 * h) % this->f);

  cout << "c1: " << c1 << endl;

  ct[0] = c0;
  ct[1] = c1;
}

void Cryptosystem::Decrypt(ZZ_pX &m, ZZ_pX ct[]) {
  ZZ_pX m_tilde;
  m_tilde.SetLength(this->n);
  m.SetLength(this->n);

  m_tilde = ct[0] + ((ct[1] * this->sk) % this->f);

  cout << "m_tilde: " << m_tilde << endl;

  cout << "m_tilde dec: ";

  for(int i = 0; i < this->n; i++) {
    ZZ x;
    conv(x, m_tilde[i]);
    if(x > to_ZZ((this->q-1)/2)) x = x - this->q;
    cout << x % 2 << " ";
  }

  cout << endl;
}

void Cryptosystem::Decrypt3(ZZ_pX &m, ZZ_pX ct[]) {
  ZZ_pX m_tilde, sksq;
  m_tilde.SetLength(this->n);
  m.SetLength(this->n);

  sksq = (this->sk * this->sk) % this->f;

  m_tilde = ct[0] + ((ct[1] * this->sk) % this->f) + ((ct[2] * sksq) % this->f);

  cout << "m_tilde: " << m_tilde << endl;

  cout << "m_tilde dec: ";

  for(int i = 0; i < this->n; i++) {
    ZZ x;
    conv(x, m_tilde[i]);
    if(x > to_ZZ((this->q-1)/2)) x = x - this->q;
    cout << x % 2 << " ";
  }

  cout << endl;
  /*ZZX r2_2, mulc1_r2_2, mulc2_r2, mm_prime, mm_d;

  PolyProd(r2_2, this->sk, this->sk);
  cout << "r2_2: " << r2_2 << endl;

  PolyProd(mulc1_r2_2, ct[0], r2_2);
  cout << "mulc1_r2_2: " << mulc1_r2_2 << endl;

  PolyProd(mulc2_r2, ct[1], this->sk);
  cout << "mulc2_r2: " << mulc2_r2 << endl;

  PolyAdd(mm_prime, mulc1_r2_2, mulc2_r2);
  cout << "mm_prime: " << mm_prime << endl;
  PolyAdd(mm_prime, mm_prime, ct[2]);
  cout << "mm_prime: " << mm_prime << endl;

  Decode(m, mm_prime);
  cout << "mm_d: " << m << endl;*/
}

Cryptosystem::Cryptosystem() {
}

void Cryptosystem::HEAdd(ZZ_pX sum[], ZZ_pX ct1[], ZZ_pX ct2[]) {

  sum[0] = ct1[0] + ct2[0];
  sum[1] = ct1[1] + ct2[1];

}

void Cryptosystem::NOT(ZZ_pX ct[]) {

  /*for(int i = 0; i < this->n; i++) {
    ct[0][i] = (ct[0][i] + ((this->q - 1)*(5/8))) % this->q;
  }

  for(int i = 0; i < this->n; i++) {
    ct[1][i] = (ct[1][i] + ((this->q - 1)*(5/8))) % this->q;
  }*/

}

void Cryptosystem::AND(ZZ_pX result[], ZZ_pX ct1[], ZZ_pX ct2[]) {
  /*result[0].SetLength(this->n);
  result[1].SetLength(this->n);

  for(int i = 0; i < this->n; i++) {
    result[0][i] = (ct1[0][i] & ct2[0][i]);// % this->q;
    result[1][i] = (ct1[1][i] & ct2[1][i]);// % this->q;
  }*/

}

void Cryptosystem::HEMult(ZZ_pX prod[], ZZ_pX ct1[], ZZ_pX ct2[]) {
  prod[0] = (ct1[0] * ct2[0]) % this->f;
  prod[1] = ((ct1[0] * ct2[1]) % this->f) + ((ct1[1] * ct2[0]) % this->f);
  prod[2] = (ct1[1] * ct2[1]) % this->f;
}

void Cryptosystem::PolyRandom(ZZX& poly, long mod) {
  /*poly.SetLength(this->n);
  std::random_device rd;
  std::mt19937 gen{rd()};
  std::uniform_int_distribution<> dis{0, INT_MAX};

  for(int i = 0; i < this->n; i++) {
    long coeff = dis(gen);
    //std::cout << "coeff: " << coeff <<std::endl;
    poly[i] = (coeff) % mod;
  }*/
}

void Cryptosystem::PolyRandom(ZZX& poly) {
  /*PolyRandom(poly, 7);*/
}

void Cryptosystem::PolyAdd(ZZX& sum, ZZX& a, ZZX& b) {
  /*sum.SetLength(this->n);
  for(int i = 0; i < this->n; i++)
    sum[i] = (a[i] + b[i]) % this->q;*/
}

void Cryptosystem::PolySub(ZZX& sub, ZZX& a, ZZX& b) {
  /*sub.SetLength(this->n);
  for(int i = 0; i < this->n; i++)
    sub[i] = (a[i] - b[i]) % this->q;*/
}

void Cryptosystem::PolyProdByScalar(ZZX& prod, ZZX& a, ZZ b) {
  /*prod.SetLength(this->n);

  prod = a * b;

  for(int i = 0; i < this->n; i++)
    prod[i] = (prod[i] % this->q) % this->q;

  prod = prod % this->f;*/
}

void Cryptosystem::PolyProd(ZZX& prod, ZZX& a, ZZX& b) {
  /*prod.SetLength(this->n);
  NTL::MulMod(prod, a, b, this->f);

  for(int i = 0; i < this->n; i++)
    prod[i] = prod[i] % this->q;*/
}

void Cryptosystem::Encode(ZZX& m_bar, ZZX& m) {

  /*m_bar.SetLength(this->n);

  for(int i = 0; i < this->n; i++)
    m_bar[i] = m[i] > 0 ? to_ZZ((long)((this->q - 1) / 2)) : to_ZZ(0);*/
}

void Cryptosystem::Decode(ZZX& m_d, ZZX& m_prime) {
  /*m_d.SetLength(this->n);
  ZZ lowerb = to_ZZ((this->q-1)/4);
  ZZ upperb = to_ZZ(3*lowerb);

  for(int i = 0; i < this->n; i++) {
    cout << "m_prime[i]: " << m_prime[i] << ", lowerb: " << lowerb << ", upperb: " << upperb << endl;
      if(m_prime[i] >= lowerb && m_prime[i] < upperb)
          m_d[i] = (int32_t) 1;
      else
          m_d[i] = (int32_t) 0;
  }*/
}
