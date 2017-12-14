/*int k = 3, q = 97, t = 2;
int n = 8; //n = 2^k
ZZX s, f, e;

f.SetLength(n + 1);
f[8] = 1;
f[0] = 1;

s.SetLength(n);
//s[7] = -1; s[6] = -1; s[5] = -1; s[4] = 1; s[3] = 1; s[2] = 1; s[1] = 1; s[0] = -1;
s[7] = 96; s[6] = 96; s[5] = 96; s[4] = 1; s[3] = 1; s[2] = 1; s[1] = 1; s[0] = 96;

e.SetLength(n);
//e[7] = 0; e[6] = -2; e[5] = 0; e[4] = 0; e[3] = -2; e[2] = 2; e[1] = -1; e[0] = 1;
e[7] = 0; e[6] = 95; e[5] = 0; e[4] = 0; e[3] = 95; e[2] = 2; e[1] = 96; e[0] = 1;
//e[7] = 1; e[6] = 2; e[5] = 1; e[4] = 2; e[3] = 1; e[2] = 2; e[1] = 1; e[0] = 1;

ZZX a0, a1;

a1.SetLength(n);
//a1[7] = 27; a1[6] = -11; a1[5] = -33; a1[4] = 41; a1[3] = -18; a1[2] = -5; a1[1] = -37; a1[0] = -16;
a1[7] = 27; a1[6] = 86; a1[5] = 64; a1[4] = 41; a1[3] = 79; a1[2] = 92; a1[1] = 60; a1[0] = 81;

a0.SetLength(n);
//a0[7] = 10; a0[6] = -25; a0[5] = 46; a0[4] = -37; a0[3] = 23; a0[2] = 27; a0[1] = -43; a0[0] = 31;
a0[7] = 10; a0[6] = 72; a0[5] = 46; a0[4] = 60; a0[3] = 23; a0[2] = 27; a0[1] = 54; a0[0] = 31;

ZZX pk[2], pk0_0, pk0_1;
Prod(pk0_0, a1, s, f, to_ZZ(q));
ProdByScalar(pk0_1, e, to_ZZ(t), f, to_ZZ(q));
Add(pk[0], pk0_0, pk0_1, to_ZZ(q));
pk[1] = a1;

cout << "pk_0: " << pk[0] << endl;
cout << "pk_1: " << pk[1] << endl;

ZZX m, ct[2], u, g, h, ct0_1, ct0_2, ct_partial, ct1_1, ct1_2;
m.SetLength(n);
m[7] = 0; m[6] = 1; m[5] = 0; m[4] = 1; m[3] = 1; m[2] = 0; m[1] = 0; m[0] = 1;

u.SetLength(n);
u[7] = 0; u[6] = -2; u[5] = 3; u[4] = 0; u[3] = 2; u[2] = 0; u[1] = -1; u[0] = 0;

g.SetLength(n);
g[7] = 0; g[6] = -1; g[5] = 0; g[4] = 0; g[3] = 0; g[2] = -1; g[1] = 2; g[0] = 0;

h.SetLength(n);
h[7] = -1; h[6] = 0; h[5] = 1; h[4] = 1; h[3] = 0; h[2] = 0; h[1] = 1; h[0] = 1;

Prod(ct0_1, pk[0], u, f, to_ZZ(q));
ProdByScalar(ct0_2, g, to_ZZ(2), f, to_ZZ(q));
Add(ct_partial, ct0_1, ct0_2, to_ZZ(q));
Add(ct[0], ct_partial, m, to_ZZ(q));

cout << "ct_0: " << ct[0] << endl;

Prod(ct1_1, pk[1], u, f, to_ZZ(q));
ProdByScalar(ct1_2, h, to_ZZ(2), f, to_ZZ(q));
Add(ct[1], ct1_1, ct1_2, to_ZZ(q));

cout << "ct_1: " << ct[1] << endl;

ZZX m_tilde, ct1_s;

Prod(ct1_s, ct[1], s, f, to_ZZ(q));
Add(m_tilde, ct[0], ct1_s, to_ZZ(q));

cout << "m_tilde: " << m_tilde << endl;

cout << "m: " << m << endl;

for(int i = 0; i < N; i++)
  cout << m_tilde[i] % 2 << ", ";*/
