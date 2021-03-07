/*

Copyright 2020 by Suyash Bagad, Saravanan Vijayakumaran

This file is part of MProvePlus-Ristretto library
Link: https://github.com/suyash67/MProvePlus-Ristretto

*/

#![allow(non_snake_case)]

use alloc::string::String;
use Errors::{self, OmniresError};
use crate::errors::ProofError;

use curve25519_dalek::ristretto::{RistrettoPoint};
use curve25519_dalek::traits::VartimeMultiscalarMul;
use curve25519_dalek::scalar::Scalar;
use crate::proofs::inner_product::InnerProductProof;
use crate::util;
use alloc::vec::Vec;
use alloc::vec;
use merlin::Transcript;
use core::iter;
use sha2::Sha512;
use core::cmp;
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use curve25519_dalek::constants;
use crate::generators::{PedersenGens};
use alloc::borrow::Borrow;
use libc_print::{libc_println};
use to_binary::{BinaryString,BinaryError};


static NUM_BITS: usize = 64;

#[derive(Clone, Debug)]
pub struct Constraints{
    alpha: Vec<Scalar>,
    beta: Vec<Scalar>,
    theta: Vec<Scalar>,
    theta_inv: Vec<Scalar>,
    nu: Vec<Scalar>,
    zeta: Vec<Scalar>,
    delta: Scalar,
}

impl Constraints{
    pub fn generate_constraints(
        u: Scalar,
        v: Scalar,
        y: Scalar,
        z: Scalar,
        n: usize,
        s: usize,
        // left: Vec<Scalar>,
        // right: Vec<Scalar>
    ) -> Constraints {

        // vector sizes
        let t: usize = s*n + n + 3*s + NUM_BITS + 2;
        let sn: usize = s*n;

        let one = Scalar::one();
        let two = Scalar::from(2 as u32);

        let y_s: Vec<Scalar> = util::exp_iter(y).take(s).collect();
        let y_n: Vec<Scalar> = util::exp_iter(y).take(n).collect();
        let y_sn_plus_bits: Vec<Scalar> = util::exp_iter(y).take(sn + NUM_BITS).collect();
        let v_s: Vec<Scalar> = util::exp_iter(v).take(s).collect();
        let two_bits: Vec<Scalar> = util::exp_iter(two).take(NUM_BITS).collect();
        
        // v^s ⊗ y^n 
        let mut vs_kronecker_yn: Vec<Scalar> = Vec::with_capacity(sn);
        vs_kronecker_yn.extend_from_slice(&y_n.as_slice());
        
        // 1^s ⊗ y^n 
        let mut ones_kronecker_yn: Vec<Scalar> = Vec::with_capacity(sn);
        ones_kronecker_yn.extend_from_slice(y_n.as_slice());

        // y^s ⊗ 1^n 
        let mut ys_kronecker_one: Vec<Scalar> = Vec::with_capacity(sn);
        ys_kronecker_one.extend_from_slice(&vec![Scalar::one(); n]);
        let mut temp_vec;
        for i in 1..s {
            temp_vec = y_n.iter().map(|yi| yi * v_s[i]).collect::<Vec<Scalar>>();
            vs_kronecker_yn.extend_from_slice(&temp_vec);
            ones_kronecker_yn.extend_from_slice(&y_n);
            ys_kronecker_one.extend_from_slice(&vec![y_s[i]; n]);
        }

        // v0 = [1  1  1^n  0^s  y^{sn + bits}  0^{2s}]
        let mut v0: Vec<Scalar> = Vec::with_capacity(t);
        v0.extend_from_slice(&vec![one; n+2]);
        v0.extend_from_slice(&vec![Scalar::zero(); s]);
        v0.extend_from_slice(&y_sn_plus_bits);
        v0.extend_from_slice(&vec![one; 2*s]);

        // v1 = [0  0  0^n  y^s  0^{sn + bits + 2s}]
        let mut v1: Vec<Scalar> = Vec::with_capacity(t);
        v1.extend_from_slice(&vec![Scalar::zero(); n+2]);
        v1.extend_from_slice(&y_s);
        v1.extend_from_slice(&vec![Scalar::zero(); sn + NUM_BITS + 2*s]);

        // v2 = [0  0  0^n  0^s  0^{sn}  2^{bits}  0^{2s}]
        let mut v2: Vec<Scalar> = Vec::with_capacity(t);
        v2.extend_from_slice(&vec![Scalar::zero(); sn+n+s+2]);
        v2.extend_from_slice(&two_bits);
        v2.extend_from_slice(&vec![Scalar::zero(); 2*s]);

        // v3 = [0  0  0^n  0^s  (y^s ⊗ 1^n)  0^{bits}  0^{2s}]
        let mut v3: Vec<Scalar> = Vec::with_capacity(t);
        v3.extend_from_slice(&vec![Scalar::zero(); n+s+2]);
        v3.extend_from_slice(&ys_kronecker_one);
        v3.extend_from_slice(&vec![Scalar::zero(); NUM_BITS + 2*s]);

        // v4 = [1  0  0^n  0^s  0^{sn}  0^{bits}  u.v^s  0^s]
        let mut v4: Vec<Scalar> = Vec::with_capacity(t);
        let u_v_s: Vec<Scalar> = v_s.iter().zip(vec![u; s].iter()).map(|(vi, u)| vi * u).collect();
        v4.push(Scalar::one());
        v4.extend_from_slice(&vec![Scalar::zero(); sn+n+s+1 + NUM_BITS]);
        v4.extend_from_slice(&u_v_s);
        v4.extend_from_slice(&vec![Scalar::zero(); s]);

        // v5 = [0  1  0^n  0^s  0^{sn}  0^{bits}  0^s  u.v^s]
        let mut v5: Vec<Scalar> = Vec::with_capacity(t);
        v5.push(Scalar::zero());
        v5.push(Scalar::one());
        v5.extend_from_slice(&vec![Scalar::zero(); sn+n+ 2*s + NUM_BITS]);
        v5.extend_from_slice(&u_v_s);

        // v6 = [0  0  -y^n  0^s  (v^s ⊗ y^n)  0^{bits}  0^s  0^s]
        let mut v6: Vec<Scalar> = Vec::with_capacity(t);
        let minus_y_n: Vec<Scalar> = (0..n).map(|i| -y_n[i]).collect();
        v6.extend_from_slice(&vec![Scalar::zero(); 2]);
        v6.extend_from_slice(&minus_y_n);
        v6.extend_from_slice(&vec![Scalar::zero(); s]);
        v6.extend_from_slice(&vs_kronecker_yn);
        v6.extend_from_slice(&vec![Scalar::zero(); NUM_BITS + 2*s]);

        // v7 = [0  0  0^n  0^s  0^{sn}  2^{bits}  -1^s  0^s]
        let mut v7: Vec<Scalar> = Vec::with_capacity(t);
        v7.extend_from_slice(&vec![Scalar::zero(); sn+n+s+2]);
        v7.extend_from_slice(&two_bits);
        v7.extend_from_slice(&vec![-Scalar::one(); s]);
        v7.extend_from_slice(&vec![Scalar::zero(); s]);

        // v8 = [0  0  0^n  0^s  y^{sn + bits}  0^{bits}  0^s  0^s]
        let mut v8: Vec<Scalar> = Vec::with_capacity(t);
        v8.extend_from_slice(&vec![Scalar::zero(); n+s+2]);
        v8.extend_from_slice(&y_sn_plus_bits);
        v8.extend_from_slice(&vec![Scalar::zero(); NUM_BITS + 2*s]);

        // u5 = [0  0  0^n  v^s  0^{sn}  0^{bits}  0^s  0^s]
        let mut u5: Vec<Scalar> = Vec::with_capacity(t);
        u5.extend_from_slice(&vec![Scalar::zero(); n+2]);
        u5.extend_from_slice(&v_s);
        u5.extend_from_slice(&vec![Scalar::zero(); sn + 2*s + NUM_BITS]);

        // theta = v0 + z.v1
        let theta: Vec<Scalar> = v0.iter().zip(v1.iter()).map(|(v0i, v1i)| v0i + z * v1i).collect();

        // theta_inv
        let mut theta_inv = theta.clone();
        util::batch_invert(&mut theta_inv);

        // zeta = z^2.v2 + z^3.v3 + ... + z^7.v7
        let z2 = z * z;
        let z3 = z * z2;
        let z4 = z2 * z2;
        let z5 = z * z4;
        let z6 = z2 * z4;
        let z7 = z * z6;
        let z8 = z * z7;
        let zeta: Vec<Scalar> = (0..t).map(|i| 
                z2*v2[i] + z3*v3[i] + z4*v4[i] + z5*v5[i] + z6*v6[i] + z7*v7[i]
            ).collect();

        // nu = z8.v8
        let nu: Vec<Scalar> = (0..t).map(|i| z8*v8[i]).collect();

        // alpha = theta^{-1} ∘ (nu + z5.u5)
        // No (-) sign before nu since we want <cL+cR-1^{t}, v6>=0
        let alpha: Vec<Scalar> = (0..t).map(|i| 
                theta_inv[i] * (z5*u5[i] + nu[i])
            ).collect();

        // beta = theta^{-1} ∘ (zeta + nu)
        let beta: Vec<Scalar> = (0..t).map(|i| 
                theta_inv[i] * (zeta[i] + nu[i])
            ).collect();

        // delta = (z + z3).(1 + y + ... + y^{s-1}) + <1^t, nu> + <alpha, zeta + nu>
        let ones_ys: Scalar = y_s.iter().sum();
        let onet_nu: Scalar = nu.iter().sum();
        let delta_cons: Scalar = (z + z3)*ones_ys + onet_nu;
        let delta: Scalar = (0..t).map(|i| 
                alpha[i] * (zeta[i] + nu[i])
            ).fold(delta_cons, |acc, x| acc + x);
        
        return Constraints{
            alpha,
            beta,
            theta,
            theta_inv,
            nu,
            zeta,
            delta,
        };
        
    }
}

#[derive(Clone, Debug)]
pub struct Omnires {
    // vector of key-images
    I_vec: Vec<RistrettoPoint>,
    // commitment to the total reserves
    C_res: RistrettoPoint,
    // commitment to the secret vectors
    A: RistrettoPoint,
    // commitment to the blinding factors
    S: RistrettoPoint,
    // commitment to the \\(t_1\\) coefficient of \\( t(x) \\)
    T1: RistrettoPoint,
    // commitment to the \\(t_2\\) coefficient of \\( t(x) \\)
    T2: RistrettoPoint,
    // evaluation of the polynomial \\(t(x)\\) at the challenge point \\(x\\)
    tau_x: Scalar,
    // blinding factor for the synthetic commitment to \\(t(x)\\)
    r: Scalar,
    // blinding factor for the synthetic commitment to the inner-product argument
    t_hat: Scalar,
    // proof data for the inner-product argument.
    inner_product_proof: InnerProductProof,
    // constraint vectors
    constraint_vec: Constraints,
}

impl Omnires {
    pub fn prove(
        // crs
        G: &RistrettoPoint,                  // base for amount in a pedersen commitment
        H: &RistrettoPoint,                  // base for blinding factor in a pedersen commitment
        H_prime: &RistrettoPoint,            // base used in commitment to inner product, F in the paper
        p_vec: &[RistrettoPoint],            // P vector in the paper 
        g_prime_vec: &[RistrettoPoint],      // G' vector in the paper
        h_vec: &[RistrettoPoint],            // H vector in the paper
        g_vec_append: &[RistrettoPoint],     // needed to bring G_w to a power of two
        h_vec_append: &[RistrettoPoint],     // needed to bring H to a power of two
        // stmt
        C_vec: &[RistrettoPoint],            // vector of commitments
        P_vec: &[RistrettoPoint],            // addresses in the ring (public keys), R vector in the paper
        H_vec: &[RistrettoPoint],            // hash of addresses
        // wit
        E_vec: &[Scalar],                    // secret indices
        x_vec: &[Scalar],                    // secret keys
        amounts: &[Scalar],                  // source amounts
        blindings: &[Scalar],                // blinding factors of source amounts
        a_res: &Scalar,                      // amount of total reserves 
        gamma: &Scalar,                      // blinding factor of total reserves
    ) -> Omnires {

        // Ring size and number of addresses owned by the exchange
        let n = E_vec.len();
        let s = x_vec.len();
        let sn = s*n;

        // size of honestly encoded witness vector
        let t = sn + n + 3*s + NUM_BITS + 2;
        assert_eq!(t, h_vec.len());

        // transcript used for generating challenges
        let mut transcript: Vec<u8> = Vec::with_capacity(1000);

        // other prelims
        let p_len = n + s + 2;
        let mut rng = rand::thread_rng();

        // generate u,v
        transcript.extend_from_slice(G.compress().as_bytes());
        let u = Scalar::hash_from_bytes::<Sha512>(&transcript);
        transcript.extend_from_slice(H.compress().as_bytes());
        let v = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // define E_mat and e_hat
        let mut E_mat: Vec<Scalar> = vec![Scalar::zero(); sn];
        
        let mut I_vec_temp: Vec<RistrettoPoint> = vec![*G; s];

        let v_s: Vec<Scalar> = util::exp_iter(v).take(s).collect();
        
        let mut index: usize = 0;
        let e_hat: Vec<Scalar> = (0..n)
            .map(|i| {
                if E_vec[i] == Scalar::one() {
                    // update 'index' numbered row of E_mat
                    E_mat[n*index + i] = Scalar::one();
                    
                    // generate key-images I = (x_i)H_i
                    I_vec_temp[index] = H_vec[i] * x_vec[index];

                    index = index + 1;

                    // i-th power of v
                    v_s[index-1]
                }
                else {
                    Scalar::zero()
                }
            })
            .collect();

        // E_mat_comp is !(E_mat)
        let E_mat_comp: Vec<Scalar> = (0..sn)
            .map(|i| 
                if E_mat[i]==Scalar::zero() {
                    Scalar::one()
                }
                else {
                    Scalar::zero()
                }
            )
            .collect();

        // generate B (the binary representation of a_res)
        let a_res_bytes = a_res.to_bytes();
        let non_zero_bytes = NUM_BITS / 8;
        for i in non_zero_bytes..32 {
            assert!(a_res_bytes[i] == 0u8, "Amount of total reserves is greater than 2^64!");
        }
        let a_res_binary = &BinaryString::from(a_res_bytes).0[0..NUM_BITS];

        let mut a_res_final: String = String::new();
        for i in 0..NUM_BITS/8 {
            let reversed = &a_res_binary[8*i..8*(i+1)].chars().rev().collect::<String>();
            a_res_final.push_str(reversed);
        }

        let mut B: Vec<Scalar> = Vec::with_capacity(NUM_BITS);
        let mut B_comp: Vec<Scalar> = Vec::with_capacity(NUM_BITS);
        for c in a_res_final.chars() {
            if c == '1' {
                B.push(Scalar::one());
                B_comp.push(Scalar::zero());
            }
            else {
                B.push(Scalar::zero());
                B_comp.push(Scalar::one());
            }
        }

        // generate key-images from I_vec_temp
        let I_vec = I_vec_temp;

        // generate commitment C_res to total assets
        let C_res = H * a_res + G * gamma;

        // define compressed stmt Y_hat_vec
        let u_square = u*u;
        let Y_hat_vec: Vec<RistrettoPoint> = (0..n).map(|i| 
                P_vec[i] + (C_vec[i] * u) + (H_vec[i] * u_square)
            ).collect();
        
        // define compressed I_hat_vec
        let minus_u_sq = -u_square;
        let I_hat_vec: Vec<RistrettoPoint> = (0..s).map(|i| I_vec[i] * (minus_u_sq * v_s[i])).collect();

        // define x_inv_vec
        let mut x_inv_vec = x_vec.clone().to_vec();
        util::batch_invert(&mut x_inv_vec);

        // define compressed secrets
        // xi = - <v^s, u.amounts>
        // eta = - <v^s, (x_vec + u.blindings)>
        let minus_u = -u;
        let xi = (0..s).map(|i| v_s[i] * amounts[i] * minus_u).fold(Scalar::zero(), |acc, x| acc + x);
        let eta = (0..s).map(|i| -v_s[i] * (x_vec[i] + u * blindings[i])).fold(Scalar::zero(), |acc, x| acc + x);

        // DEBUG /////////////////
        // CHECK main equality
        let main_eq = RistrettoPoint::vartime_multiscalar_mul(
            e_hat.iter()
            .chain(x_inv_vec.iter())
            .chain(iter::once(&xi))
            .chain(iter::once(&eta)),
            Y_hat_vec.iter()
            .chain(I_hat_vec.iter())
            .chain(iter::once(H))
            .chain(iter::once(G)),
        );
        assert_eq!(main_eq, Scalar::zero() * G);
        //////////////////////////

        // secret vectors
        // c_L := ( ξ || η || ê || x^◦−1 || vec(E) || vec(B) || amounts || blindings)
        let mut c_L: Vec<Scalar> = Vec::with_capacity(t);
        c_L.push(xi.clone());
        c_L.push(eta.clone());
        c_L.extend_from_slice(&e_hat.clone());
        c_L.extend_from_slice(&x_inv_vec.clone());
        c_L.extend_from_slice(&E_mat.clone());
        c_L.extend_from_slice(&B.clone());
        c_L.extend_from_slice(&amounts.clone());
        c_L.extend_from_slice(&blindings.clone());
        
        // c R := (0^{2+n} || x || vec(E) − 1^sn || vec(B) - 1^s || 0^{2s} )
        let mut c_R: Vec<Scalar> = Vec::with_capacity(t);
        c_R.extend_from_slice(&vec![Scalar::zero(); n + 2]);
        c_R.extend_from_slice(&x_vec);
        c_R.extend_from_slice(&E_mat_comp);
        c_R.extend_from_slice(&B_comp);
        c_R.extend_from_slice(&vec![Scalar::zero(); 2*s]);
    
        // defining g_0
        let mut g_vec_0: Vec<RistrettoPoint> = Vec::with_capacity(t);
        g_vec_0.extend_from_slice(&p_vec);
        g_vec_0.extend_from_slice(&g_prime_vec);

        // P -> V: A
        // A = (H' * r_a) + <g_0 * c_L> + <h * c_R>
        let r_A = Scalar::random(&mut rng);  
        let A = RistrettoPoint::vartime_multiscalar_mul(
            c_L.iter()
            .chain(c_R.iter())
            .chain(iter::once(&r_A)),
            g_vec_0.iter()
            .chain(h_vec.iter())
            .chain(iter::once(H_prime)),
        );

        // challenge w
        transcript.extend_from_slice(A.compress().as_bytes());
        let w = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // defining g_w
        let mut g_vec_w: Vec<RistrettoPoint> = Vec::with_capacity(t);
        g_vec_w.extend_from_slice(&vec![*H, *G]);
        g_vec_w.extend_from_slice(&Y_hat_vec);
        g_vec_w.extend_from_slice(&I_hat_vec);

        g_vec_w = (0..p_len)
            .map(|i| {
                p_vec[i] + (g_vec_w[i] * w)
            })
            .collect();
        g_vec_w.extend_from_slice(&g_prime_vec);
        
        // P -> V: S
        // S = (H' * r_s) + <g_w * s_L> + <h * s_R>
        let r_S = Scalar::random(&mut rng);
        let s_L = (0..t).map(|_| Scalar::random(&mut rng)).collect::<Vec<Scalar>>();

        // s_R[j]=0 if c_R[j]=0, else a random scalar
        let mut s_R: Vec<Scalar> = Vec::with_capacity(t);
        s_R.extend_from_slice(&vec![Scalar::zero(); n + 2]);
        let s_R_ext = (0..s+sn+NUM_BITS).map(|_| Scalar::random(&mut rng)).collect::<Vec<Scalar>>();
        s_R.extend_from_slice(&s_R_ext);
        s_R.extend_from_slice(&vec![Scalar::zero(); 2*s]);

        let r_start = n + 2;
        let r_end = t - 2 * s;
        let S = RistrettoPoint::vartime_multiscalar_mul(
            s_L.iter()
            .chain(s_R[r_start..r_end].iter())
            .chain(iter::once(&r_S)),
            g_vec_w.iter()
            .chain(h_vec[r_start..r_end].iter())
            .chain(iter::once(H_prime)),
        );

        // challenges y,z
        transcript.extend_from_slice(S.compress().as_bytes());
        let y = Scalar::hash_from_bytes::<Sha512>(&transcript);
        transcript.extend_from_slice(S.compress().as_bytes());
        let z = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // generate constraint vectors
        let constraint_vec = Constraints::generate_constraints(u,v,y,z,n,s);

        let theta = constraint_vec.theta.clone();
        let theta_inv = constraint_vec.theta_inv.clone();
        let alpha = constraint_vec.alpha.clone();
        let zeta = constraint_vec.zeta.clone();
        let nu = constraint_vec.nu.clone();
        
        // calculate t2, t1, t0
        let t2 = (0..t).map(|i| (s_L[i] * s_R[i]) * theta[i]).fold(Scalar::zero(), |acc, x| acc + x);
        
        let t1 = (0..t)
            .map(|i| {
                (c_L[i] + alpha[i])*theta[i]*s_R[i] + (c_R[i]*theta[i] + zeta[i] + nu[i])*s_L[i]
            })
            .fold(Scalar::zero(), |acc, x| acc + x); 

        // P -> V: T_1, T_2
        let tau1 = Scalar::random(&mut rng);
        let tau2 = Scalar::random(&mut rng);
        let T1 = H*t1 + G*tau1;
        let T2 = H*t2 + G*tau2;

        // generate challenge x
        transcript.extend_from_slice(T1.compress().as_bytes());
        transcript.extend_from_slice(T2.compress().as_bytes());
        let x = Scalar::hash_from_bytes::<Sha512>(&transcript);
        
        // P -> V: tau_x, r, t_hat
        // compute tau_x, r, Lp, Rp, t_hat
        let tau_x = z*z*gamma + x*tau1 + x*x*tau2;

        let r = (r_S*x) + r_A;
        
        let Lp: Vec<Scalar> = (0..t).map(|i| s_L[i]*x + (c_L[i] + alpha[i])).collect();

        let Rp: Vec<Scalar> = (0..t).map(|i| theta[i]*(c_R[i] + s_R[i]*x) + zeta[i] + nu[i]).collect();

        let t_hat = Lp.iter().zip(Rp.iter()).fold(Scalar::zero(), |acc, x| {
            acc + (x.0*x.1)
        });

        // DEBUG ///////////////// 
        let lhs = t2*x*x + t1*x + constraint_vec.delta + z*z*a_res;
        assert_eq!(lhs, t_hat);
        ///////////////////

        // Running inner product argument
        let Q = RistrettoPoint::hash_from_bytes::<Sha512>(b"test point");

        // Run ipp for non-power of two secret vectors
        let N = t.next_power_of_two();
        let res = N-t;
        let zero_append_vec = vec![Scalar::zero();res];

        // Append 0s to secret vectors
        let mut a: Vec<Scalar> = Vec::with_capacity(N);
        let mut b: Vec<Scalar> = Vec::with_capacity(N);
        // let mut a = Lp.clone();
        // let mut b = Rp.clone();
        a.extend_from_slice(&Lp);
        b.extend_from_slice(&Rp);
        a.extend_from_slice(&zero_append_vec);
        b.extend_from_slice(&zero_append_vec);

        // Multipliers of base vectors
        let G_factors: Vec<Scalar> = iter::repeat(Scalar::one()).take(N).collect();
        let mut H_factors: Vec<Scalar> = Vec::with_capacity(N);
        H_factors.extend_from_slice(&theta_inv);
        H_factors.extend_from_slice(&vec![Scalar::one(); res]);

        // Append random group generators to g_vec_w and hi_tag        
        let mut g_vec_long: Vec<RistrettoPoint> = Vec::with_capacity(N);
        let mut h_vec_long: Vec<RistrettoPoint> = Vec::with_capacity(N);
        g_vec_long.extend_from_slice(&g_vec_w);
        g_vec_long.extend_from_slice(&g_vec_append);
        h_vec_long.extend_from_slice(&h_vec);
        h_vec_long.extend_from_slice(&h_vec_append);

        let mut verifier = Transcript::new(b"innerproduct");
        let inner_product_proof = InnerProductProof::create(
            &mut verifier,
            &Q,
            &G_factors,
            &H_factors,
            g_vec_long.clone(),
            h_vec_long.clone(),
            a.clone(),
            b.clone(),
        );

        return Omnires {
            I_vec,
            C_res,
            A,
            S,
            T1,
            T2,
            tau_x,
            r,
            t_hat: t_hat,
            inner_product_proof,
            constraint_vec,
        };
    }

    pub fn verify(
        &self,
        // crs
        G: &RistrettoPoint,
        H: &RistrettoPoint,
        H_prime: &RistrettoPoint,
        p_vec: &[RistrettoPoint],
        g_prime_vec: &[RistrettoPoint],
        h_vec: &[RistrettoPoint],
        g_vec_append: &[RistrettoPoint],
        h_vec_append: &[RistrettoPoint],
        // stmt
        C_vec: &[RistrettoPoint], // vector of commitments
        P_vec: &[RistrettoPoint], // addresses in the ring (public keys)
        H_vec: &[RistrettoPoint], // hash of addresses
    ) -> Result<(), Errors> {

        // vector lengths
        let n = C_vec.len();
        let s = self.I_vec.len();
        let t = s*n + n + 3*s + NUM_BITS + 2;
        let p_len = n+s+2;
        let N = t.next_power_of_two();
        let res = N-t;

        // transcript initialization
        let mut transcript: Vec<u8> = Vec::with_capacity(1000);

        // re-generate challenges u, v, y, z
        transcript.extend_from_slice(G.compress().as_bytes());
        let u = Scalar::hash_from_bytes::<Sha512>(&transcript);
        transcript.extend_from_slice(H.compress().as_bytes());
        let v = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // re-generate challenge w and compute Q
        transcript.extend_from_slice(self.A.compress().as_bytes());
        let w = Scalar::hash_from_bytes::<Sha512>(&transcript);
        let Q = RistrettoPoint::hash_from_bytes::<Sha512>(b"test point");

        // re-generate challenge y, z
        transcript.extend_from_slice(self.S.compress().as_bytes());
        let y = Scalar::hash_from_bytes::<Sha512>(&transcript);
        transcript.extend_from_slice(self.S.compress().as_bytes());
        let z = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // re-generate challenge x
        transcript.extend_from_slice(self.T1.compress().as_bytes());
        transcript.extend_from_slice(self.T2.compress().as_bytes());
        let x = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // build constraint vectors
        let constraint_vec1 = Constraints::generate_constraints(u,v,y,z,n,s);
        let theta_inv = constraint_vec1.theta_inv.clone();
        let alpha = constraint_vec1.alpha.clone();
        let delta = constraint_vec1.delta.clone();
        let beta = constraint_vec1.beta.clone();  

        // verification equation #2
        // lhs
        let G_t_hat = H * self.t_hat;
        let Htau_x = G * self.tau_x;
        let left_side = G_t_hat + Htau_x;

        // rhs
        let z_sq = z * z;
        let Gdelta = H * delta;
        let Tx = self.T1 * x;
        let Tx_sq = self.T2 * (x*x);
        let right_side = Gdelta + Tx + Tx_sq + (self.C_res * z_sq);

        assert_eq!(left_side, right_side);

        // towards verification eqn #3
        // define compressed stmt Y_hat_vec
        let u_square = u*u;
        let Y_hat_vec: Vec<RistrettoPoint> = (0..n).map(|i| 
                P_vec[i] + (C_vec[i] * u) + (H_vec[i] * u_square)
            ).collect();

        // define compressed I_hat_vec
        let minus_u_sq = -u_square;
        let v_s: Vec<Scalar> = util::exp_iter(v).take(s).collect();
        let I_hat_vec: Vec<RistrettoPoint> = (0..s).map(|i| self.I_vec[i] * (minus_u_sq * v_s[i])).collect();

        // defining g_w
        let mut g_vec_w: Vec<RistrettoPoint> = Vec::with_capacity(t);
        g_vec_w.extend_from_slice(&vec![*H, *G]);
        g_vec_w.extend_from_slice(&Y_hat_vec);
        g_vec_w.extend_from_slice(&I_hat_vec);

        g_vec_w = (0..p_len)
            .map(|i| {
                p_vec[i] + (g_vec_w[i] * w)
            })
            .collect();
        g_vec_w.extend_from_slice(&g_prime_vec);


        // compute a commitment to l(x),r(x) efficiently
        let minus_r = -self.r;
        let P = RistrettoPoint::vartime_multiscalar_mul(
            alpha.iter()
            .chain(beta.iter())
            .chain(iter::once(&self.t_hat))
            .chain(iter::once(&x))
            .chain(iter::once(&minus_r))
            .chain(iter::once(&Scalar::one())),
            g_vec_w.iter()
            .chain(h_vec.iter())
            .chain(iter::once(&Q))
            .chain(iter::once(&self.S))
            .chain(iter::once(H_prime))
            .chain(iter::once(&self.A)),
        );

        // Multipliers of base vectors
        let G_factors: Vec<Scalar> = iter::repeat(Scalar::one()).take(N).collect();
        let mut H_factors: Vec<Scalar> = Vec::with_capacity(N);
        H_factors.extend_from_slice(&theta_inv);
        H_factors.extend_from_slice(&vec![Scalar::one(); res]);

        // Append random group generators (same as prover's) to g_vec_w and hi_tag        
        let mut g_vec_long: Vec<RistrettoPoint> = Vec::with_capacity(N);
        let mut h_vec_long: Vec<RistrettoPoint> = Vec::with_capacity(N);
        g_vec_long.extend_from_slice(&g_vec_w);
        g_vec_long.extend_from_slice(&g_vec_append);
        h_vec_long.extend_from_slice(&h_vec);
        h_vec_long.extend_from_slice(&h_vec_append);

        let mut verifier = Transcript::new(b"innerproduct");
        let verify = self.inner_product_proof.verify(N, &mut verifier, &G_factors,
            &H_factors, &P, &Q, &g_vec_long, &h_vec_long);

        // check all three conditions are true
        if verify.is_ok() && left_side==right_side{
            Ok(())
        } else {
            Err(OmniresError)
        }
    }

    pub fn fast_verify(
        &self,
        // crs
        G: &RistrettoPoint,
        H: &RistrettoPoint,
        H_prime: &RistrettoPoint,
        p_vec: &[RistrettoPoint],
        g_prime_vec: &[RistrettoPoint],
        h_vec: &[RistrettoPoint],
        g_vec_append: &[RistrettoPoint],
        h_vec_append: &[RistrettoPoint],
        // stmt
        C_vec: &[RistrettoPoint], // vector of commitments
        P_vec: &[RistrettoPoint], // addresses in the ring (public keys)
        H_vec: &[RistrettoPoint], // hash of addresses
    ) -> Result<(), Errors> {

        // vector lengths
        let n = C_vec.len();
        let s = self.I_vec.len();
        let t = s*n + n + 3*s + NUM_BITS + 2;
        let N = t.next_power_of_two();
        let res = N-t;

        // transcript initialization
        let mut transcript: Vec<u8> = Vec::with_capacity(1000);

        // re-generate challenges u, v, y, z
        transcript.extend_from_slice(G.compress().as_bytes());
        let u = Scalar::hash_from_bytes::<Sha512>(&transcript);
        transcript.extend_from_slice(H.compress().as_bytes());
        let v = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // re-generate challenge w and compute Q
        transcript.extend_from_slice(self.A.compress().as_bytes());
        let w = Scalar::hash_from_bytes::<Sha512>(&transcript);
        let Q = RistrettoPoint::hash_from_bytes::<Sha512>(b"test point");

        // re-generate challenge y, z
        transcript.extend_from_slice(self.S.compress().as_bytes());
        let y = Scalar::hash_from_bytes::<Sha512>(&transcript);
        transcript.extend_from_slice(self.S.compress().as_bytes());
        let z = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // re-generate challenge x
        transcript.extend_from_slice(self.T1.compress().as_bytes());
        transcript.extend_from_slice(self.T2.compress().as_bytes());
        let x = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // generate scalar c for combining verification equations
        transcript.extend_from_slice(self.T2.compress().as_bytes());
        let c = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // build constraint vectors
        let constraint_vec1 = Constraints::generate_constraints(u,v,y,z,n,s);
        let theta_inv = constraint_vec1.theta_inv.clone();
        let alpha = constraint_vec1.alpha.clone();
        let delta = constraint_vec1.delta.clone();
        let beta = constraint_vec1.beta.clone();

        // compute sg_vec and sh_vec
        let mut verifier = Transcript::new(b"innerproduct");
        let (u_sq, u_inv_sq, s_vec) = InnerProductProof::verification_scalars(
            &self.inner_product_proof, 
            N, 
            &mut verifier)
            .expect(
                "Issues in the input inner product argument encountered!"
            );

        // compute exponent of g_vec_w
        let mut alpha_ext = Vec::with_capacity(N);
        alpha_ext.extend_from_slice(&alpha);
        alpha_ext.extend_from_slice(&vec![Scalar::zero(); res]);        

        let a_s_minus_alpha: Vec<Scalar> = alpha_ext
            .into_iter()
            .zip(s_vec.iter())
            .map(|(alpha_i, s_i)| (self.inner_product_proof.a * s_i) - alpha_i.borrow())
            .take(N)
            .collect();

        let mut scalar_gw = Vec::with_capacity(N);
        scalar_gw.push(&w * a_s_minus_alpha[0]);
        scalar_gw.push(&w * a_s_minus_alpha[1]);
        
        let wu = &w * &u;
        let wu_sq = &wu * &u;
        let minus_wu_sq = -&wu_sq;
        let scalar_C_vec: Vec<Scalar> = (0..n).map(|i| &wu * a_s_minus_alpha[2+i]).collect();
        let scalar_H_vec: Vec<Scalar> = (0..n).map(|i| &wu_sq * a_s_minus_alpha[2+i]).collect();
        let scalar_P_vec: Vec<Scalar> = (0..n).map(|i| &w * a_s_minus_alpha[2+i]).collect();

        let v_s: Vec<Scalar> = util::exp_iter(v).take(s).collect();
        let scalar_I_vec: Vec<Scalar> = (0..s)
            .map(|i| {
                let minus_wu_sq_vi = &minus_wu_sq * v_s[i];
                &minus_wu_sq_vi * a_s_minus_alpha[n+2+i]
            })
            .collect();

        let scalar_p_vec: Vec<Scalar> = a_s_minus_alpha[0..(n+s+2)].to_vec();
        let scalar_g_prime_vec: Vec<Scalar> = a_s_minus_alpha[(n+s+2)..t].to_vec();
        let scalar_g_extend_vec: Vec<Scalar> = a_s_minus_alpha[t..N].to_vec();

        scalar_gw.extend_from_slice(&scalar_P_vec);
        scalar_gw.extend_from_slice(&scalar_C_vec);
        scalar_gw.extend_from_slice(&scalar_H_vec);
        scalar_gw.extend_from_slice(&scalar_I_vec);
        scalar_gw.extend_from_slice(&scalar_p_vec);
        scalar_gw.extend_from_slice(&scalar_g_prime_vec);
        scalar_gw.extend_from_slice(&scalar_g_extend_vec);

        let mut points_gw: Vec<RistrettoPoint> = Vec::with_capacity(N);
        points_gw.extend_from_slice(&vec![*H, *G]);
        points_gw.extend_from_slice(&P_vec);
        points_gw.extend_from_slice(&C_vec);
        points_gw.extend_from_slice(&H_vec);
        points_gw.extend_from_slice(&self.I_vec);
        points_gw.extend_from_slice(&p_vec);
        points_gw.extend_from_slice(&g_prime_vec);
        points_gw.extend_from_slice(&g_vec_append);


        // compute exponent of h_vec
        let mut beta_ext = Vec::with_capacity(N);
        beta_ext.extend_from_slice(&beta);
        beta_ext.extend_from_slice(&vec![Scalar::zero(); res]);

        // 1/s[i] is s[!i], and !i runs from n-1 to 0 as i runs from 0 to n-1
        let inv_s = s_vec.iter().rev();

        let mut H_factors: Vec<Scalar> = Vec::with_capacity(N);
        H_factors.extend_from_slice(&theta_inv);
        H_factors.extend_from_slice(&vec![Scalar::one(); res]);

        let h_times_b_div_s = H_factors
            .into_iter()
            .zip(inv_s)
            .map(|(h_i, s_i_inv)| (self.inner_product_proof.b * s_i_inv) * h_i.borrow());

        let scalar_h_vec = beta_ext
            .into_iter()
            .zip(h_times_b_div_s)
            .map(|(beta_i, h_b_s_i)| h_b_s_i - beta_i.borrow())
            .take(N);

        let mut points_h_vec: Vec<RistrettoPoint> = Vec::with_capacity(N);
        points_h_vec.extend_from_slice(&h_vec);
        points_h_vec.extend_from_slice(&h_vec_append);

        // exponents of L_vec, R_vec
        let neg_u_sq = u_sq.iter().map(|ui| -ui);
        let neg_u_inv_sq = u_inv_sq.iter().map(|ui| -ui);

        let Ls = self.inner_product_proof
            .L_vec
            .iter()
            .map(|p| p.decompress().ok_or(ProofError::VerificationError))
            .collect::<Result<Vec<_>, _>>()
            .expect("Unable to decompress L!");

        let Rs = self.inner_product_proof
            .R_vec
            .iter()
            .map(|p| p.decompress().ok_or(ProofError::VerificationError))
            .collect::<Result<Vec<_>, _>>()
            .expect("Unable to decompress R!");

        // exponent of Q, G, H_prime, H
        let scalar_Q = (self.inner_product_proof.a * self.inner_product_proof.b) - self.t_hat;
        let scalar_H = &c * (&delta - &self.t_hat);
        let scalar_H_prime = self.r;
        let scalar_G = -c * self.tau_x;

        // exponent of T_1, T_2, A, S, C_res
        let scalar_T1 = &c * &x;
        let scalar_T2 = &scalar_T1 * &x;
        // let scalar_A = -Scalar::one();
        let scalar_S = -&x;
        let scalar_C_res = &z * &z * &c;

        let expected = RistrettoPoint::vartime_multiscalar_mul(
            iter::once(scalar_S)
                .chain(scalar_gw)
                .chain(scalar_h_vec)
                .chain(neg_u_sq)
                .chain(neg_u_inv_sq)
                .chain(vec![
                    scalar_Q, 
                    scalar_G, 
                    scalar_H_prime, 
                    scalar_H, 
                    scalar_T1,
                    scalar_T2,
                    scalar_C_res,
                    ]),
            iter::once(self.S)
                .chain(points_gw)
                .chain(points_h_vec)
                .chain(Ls)
                .chain(Rs)
                .chain(vec![
                    Q,
                    *G,
                    *H_prime,
                    *H,
                    self.T1,
                    self.T2,
                    self.C_res,
                ])
        );
        
        // check the single multi-exp is equal to A
        //
        // g_w^{a.sg_vec - alpha} . G^{b.theta_inv \circ sh_vec - beta} . 
        // L_vec^{-u_sq_vec} . R_vec^{-u_inv_sq_vec} . Q^{a.b - t_hat} .
        // H^{delta - t_hat} . (H')^{r} . H^{-c tau_x} . 
        // T1^{cx} . T2^{cx^2} . S^{-x} . 
        // C_res^{cz^2} 
        // =?
        // A  
        if expected == self.A {
            Ok(())
        } else {
            Err(OmniresError)
        }
    }

    pub fn gen_params(n: usize, s: usize) -> (
        RistrettoPoint,
        RistrettoPoint,
        RistrettoPoint, 
        Vec<RistrettoPoint>,
        Vec<RistrettoPoint>,
        Vec<RistrettoPoint>,
        Vec<RistrettoPoint>,
        Vec<RistrettoPoint>,
        Vec<RistrettoPoint>,
        Vec<RistrettoPoint>,
        Vec<RistrettoPoint>,
        Vec<Scalar>,
        Vec<Scalar>,
        Vec<Scalar>,
        Vec<Scalar>,
        Scalar,
        Scalar,
    ) {
        
        let sn = s * n;
        let t = sn + n + 3*s + NUM_BITS + 2;
        let one = Scalar::one();
        let mut rng = rand::thread_rng();

        // generate random amounts in range {0,..,2^{32}-1}
        let a_vec: Vec<Scalar> = (0..s).map(|_| Scalar::from(rng.gen::<u32>())).collect();
        let a_res = a_vec.iter().sum();
        
        // generate blinding factors
        let r_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();

        // generate secret keys
        let x_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();
        
        // G, H, H' - curve points for generating outputs and key-images
        let G = constants::RISTRETTO_BASEPOINT_POINT;
        let H = RistrettoPoint::hash_from_bytes::<Sha512>(G.compress().as_bytes());
        let H_prime = RistrettoPoint::hash_from_bytes::<Sha512>(b"h_prime");

        let pgens = PedersenGens{B: H, B_blinding: G};

        // generate p_vec, g_prime_vec, h_vec
        let p_len = n+s+2;
        let g_prime_len = t-p_len;
        let h_len = t;
        let N = (t as u64).next_power_of_two();
        let res = (N as usize)-t;

        let g_prime_vec: Vec<RistrettoPoint> = (0..g_prime_len).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let h_vec: Vec<RistrettoPoint> = (0..h_len).map(|_| RistrettoPoint::random(&mut rng)).collect();

        let p_vec: Vec<RistrettoPoint> = (0..p_len).map(|_| RistrettoPoint::random(&mut rng)).collect();

        let g_vec_append: Vec<RistrettoPoint> = (0..res).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let h_vec_append: Vec<RistrettoPoint> = (0..res).map(|_| RistrettoPoint::random(&mut rng)).collect();
        
        // generate random ring addresses and Hash of them 
        let mut P_vec: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H_vec: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        
        // Select random outputs owned by the exchange
        let mut C_vec_mut: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();        
        
        // generate random index vector of size s
        let setsize = n / s;
        let mut start_idx = 0;
        let mut end_idx = cmp::max(1, setsize-1);
        let mut dist1 = Uniform::from(0..2);
        let idx = (0..s).map(|_| {
            
            dist1 = Uniform::from(start_idx..end_idx);
            start_idx = setsize + start_idx;
            end_idx =  cmp::min(n-1, end_idx + setsize);

            dist1.sample(&mut rng)
        })
        .collect::<Vec<usize>>();

        // generate commitment to total reserves
        let gamma = Scalar::random(&mut rng);
        
        let mut index = 0;
        let E_vec = (0..n)
            .map(|i| {
                if index < idx.len() {
                    if i == idx[index] {
                        // generate commitments using a_vec, r_vec
                        C_vec_mut[i as usize] = pgens.commit(a_vec[index as usize], r_vec[index as usize]);
                        P_vec[i as usize] = G * x_vec[index];
                        index = index + 1;
                        one.clone()
                    }
                    else {
                        Scalar::zero()
                    }
                }
                else{
                    Scalar::zero()
                }
            })
            .collect::<Vec<Scalar>>();

        (G, H, H_prime, 
         p_vec, 
         g_prime_vec, 
         h_vec, 
         g_vec_append, 
         h_vec_append, 
         C_vec_mut, 
         P_vec, 
         H_vec, 
         E_vec, 
         x_vec, 
         a_vec,
         r_vec,
         a_res,
         gamma)
    }
}


#[cfg(test)]
mod omniring_tests {
    
    use super::*;
    use core::cmp;
    use rand::distributions::{Distribution, Uniform};
    use rand::Rng;
    use curve25519_dalek::constants;
    use crate::generators::{PedersenGens};
    use libc_print::{libc_println};
    use time::PreciseTime;

    #[test]
    pub fn test_omnires_aok_1_in_8(){
        
        let n=8;
        let s=1;
        let sn = s * n;
        let t = sn + n + 3*s + NUM_BITS + 2;
        let one = Scalar::one();
        let mut rng = rand::thread_rng();

        // generate random amounts in range {0,..,2^{32}-1}
        let a_vec: Vec<Scalar> = (0..s).map(|_| Scalar::from(rng.gen::<u32>())).collect();
        let a_res = a_vec.iter().sum();
        
        // generate blinding factors
        let r_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();

        // generate secret keys
        let x_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();
        
        // G, H, H' - curve points for generating outputs and key-images
        let G = constants::RISTRETTO_BASEPOINT_POINT;
        let H = RistrettoPoint::hash_from_bytes::<Sha512>(G.compress().as_bytes());
        let H_prime = RistrettoPoint::hash_from_bytes::<Sha512>(b"h_prime");

        let pgens = PedersenGens{B: H, B_blinding: G};

        // generate p_vec, g_prime_vec, h_vec
        let p_len = n+s+2;
        let g_prime_len = t-p_len;
        let h_len = t;
        let N = (t as u64).next_power_of_two();
        let res = N - (t as u64);

        let g_prime_vec: Vec<RistrettoPoint> = (0..g_prime_len).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let h_vec: Vec<RistrettoPoint> = (0..h_len).map(|_| RistrettoPoint::random(&mut rng)).collect();

        let p_vec: Vec<RistrettoPoint> = (0..p_len).map(|_| RistrettoPoint::random(&mut rng)).collect();

        let g_vec_append: Vec<RistrettoPoint> = (0..res).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let h_vec_append: Vec<RistrettoPoint> = (0..res).map(|_| RistrettoPoint::random(&mut rng)).collect();
        
        // generate random ring addresses and Hash of them 
        let mut P_vec: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        let H_vec: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();
        
        // Select random outputs owned by the exchange
        let mut C_vec_mut: Vec<RistrettoPoint> = (0..n).map(|_| RistrettoPoint::random(&mut rng)).collect();        

        // generate commitment to total reserves
        let gamma = Scalar::random(&mut rng);
        
        let idx = vec![6];
        let mut index = 0;
        let E_vec = (0..n)
            .map(|i| {
                if index < idx.len() {
                    if i == idx[index] {
                        // generate commitments using a_vec, r_vec
                        C_vec_mut[i as usize] = pgens.commit(a_vec[index as usize], r_vec[index as usize]);
                        P_vec[i as usize] = G * x_vec[index];
                        index = index + 1;
                        one.clone()
                    }
                    else {
                        Scalar::zero()
                    }
                }
                else{
                    Scalar::zero()
                }
            })
            .collect::<Vec<Scalar>>();

        let mprove_test = Omnires::prove(&G, 
                                            &H, 
                                            &H_prime, 
                                            &p_vec, 
                                            &g_prime_vec, 
                                            &h_vec, 
                                            &g_vec_append, 
                                            &h_vec_append, 
                                            &C_vec_mut, 
                                            &P_vec, 
                                            &H_vec, 
                                            &E_vec, 
                                            &x_vec,
                                            &a_vec,
                                            &r_vec,
                                            &a_res, 
                                            &gamma);

        let result = mprove_test.verify(&G, 
                                        &H, 
                                        &H_prime, 
                                        &p_vec, 
                                        &g_prime_vec, 
                                        &h_vec, 
                                        &g_vec_append, 
                                        &h_vec_append, 
                                        &C_vec_mut, 
                                        &P_vec, 
                                        &H_vec);
        assert!(result.is_ok());
    }

    pub fn test_omnires_helper(n: usize, s: usize){

            let (G, H, H_prime, 
                 p_vec, 
                 g_prime_vec, 
                 h_vec, 
                 g_vec_append, 
                 h_vec_append, 
                 C_vec_mut, 
                 P_vec, 
                 H_vec, 
                 E_vec, 
                 x_vec, 
                 a_vec,
                 r_vec,
                 a_res,
                 gamma) = Omnires::gen_params(n, s);

            // libc_println!("For n={}, s={},", n, s);
            // let start = PreciseTime::now();
            let mprove_test = Omnires::prove(&G, 
                                                &H, 
                                                &H_prime, 
                                                &p_vec, 
                                                &g_prime_vec, 
                                                &h_vec, 
                                                &g_vec_append, 
                                                &h_vec_append, 
                                                &C_vec_mut, 
                                                &P_vec, 
                                                &H_vec, 
                                                &E_vec, 
                                                &x_vec, 
                                                &a_vec,
                                                &r_vec, 
                                                &a_res,
                                                &gamma);
            // let end = PreciseTime::now();
            // libc_println!("Generation time: {:?}", start.to(end));
            
            // let start = PreciseTime::now();
            // let mut transcript = Transcript::new(b"mprovetest");
            let result = mprove_test.fast_verify(&G, 
                                                 &H, 
                                                 &H_prime, 
                                                 &p_vec, 
                                                 &g_prime_vec, 
                                                 &h_vec, 
                                                 &g_vec_append, 
                                                 &h_vec_append, 
                                                 &C_vec_mut, 
                                                 &P_vec, 
                                                 &H_vec);
            // let end = PreciseTime::now();
            // libc_println!("Verification time: {:?}", start.to(end));

            assert!(result.is_ok());
    }

    #[test]
    fn omnires_2_in_100(){
        test_omnires_helper(100, 2);
    }
    
    #[test]
    fn omnires_5_in_100(){
        test_omnires_helper(100, 5);
    }

    #[test]
    fn omnires_4_in_500(){
        test_omnires_helper(500, 4);
    }

    #[test]
    fn omnires_19_in_20(){
        test_omnires_helper(20, 19);
    }

}
