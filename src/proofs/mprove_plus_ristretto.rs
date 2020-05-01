/*

Copyright 2020 by Suyash Bagad, Saravanan Vijayakumaran

This file is part of MProvePlus-Ristretto library
Link: https://github.com/suyash67/MProvePlus-Ristretto

*/

#![allow(non_snake_case)]

use Errors::{self, MProvePlusError};

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
    ) -> Constraints {

        // vector sizes
        let t: usize = s*n + 2*n + s + 3;
        let sn: usize = s*n;

        let y_s: Vec<Scalar> = util::exp_iter(y).take(s).collect();
        let y_sn: Vec<Scalar> = util::exp_iter(y).take(sn).collect();
        let y_n: Vec<Scalar> = util::exp_iter(y).take(n).collect();
        let v_s: Vec<Scalar> = util::exp_iter(v).take(s).collect();
        let mut vs_kronecker_yn: Vec<Scalar> = y_n.clone();
        let mut ones_kronecker_yn: Vec<Scalar> = y_n.clone();
        let mut ys_kronecker_one: Vec<Scalar> = vec![Scalar::one(); n];
        for i in 1..s {
            let temp_vec: Vec<Scalar> = y_n.iter().map(|yi| yi * v_s[i]).collect();
            vs_kronecker_yn.extend_from_slice(&temp_vec);
            ones_kronecker_yn.extend_from_slice(&y_n);
            ys_kronecker_one.extend_from_slice(&vec![y_s[i]; n]);
        }

        // v0
        let mut v0: Vec<Scalar> = vec![v; 2*n+3];
        v0.extend_from_slice(&vec![Scalar::zero(); s]);
        v0.extend_from_slice(&y_sn);

        // v1
        let mut v1: Vec<Scalar> = vec![Scalar::zero(); 2*n+3];
        v1.extend_from_slice(&y_s);
        v1.extend_from_slice(&vec![Scalar::zero(); sn]);

        // v2
        let mut v2: Vec<Scalar> = vec![Scalar::one(); 1];
        v2.extend_from_slice(&vec![Scalar::zero(); t-1]);

        // v3
        let minus_y_n: Vec<Scalar> = (0..n).map(|i| -y_n[i]).collect();
        let mut v3: Vec<Scalar> = vec![Scalar::zero(); 3];
        v3.extend_from_slice(&minus_y_n);
        v3.extend_from_slice(&vec![Scalar::zero(); n+s]);
        v3.extend_from_slice(&vs_kronecker_yn);

        // v4
        let mut v4: Vec<Scalar> = vec![Scalar::zero(); n+3];
        v4.extend_from_slice(&minus_y_n);
        v4.extend_from_slice(&vec![Scalar::zero(); s]);
        v4.extend_from_slice(&ones_kronecker_yn);

        // v5
        let minus_y_pow_s = -util::scalar_exp_vartime(&y, s as u64);
        let mut v5: Vec<Scalar> = vec![Scalar::zero(), minus_y_pow_s];
        v5.extend_from_slice(&vec![Scalar::zero(); 2*n+s+1]);
        v5.extend_from_slice(&ys_kronecker_one);

        // v6
        let mut v6: Vec<Scalar> = vec![Scalar::zero(); 2*n + s + 3];
        v6.extend_from_slice(&y_sn);

        // v7
        let u_v_s: Vec<Scalar> = v_s.iter().zip(vec![u; s].iter()).map(|(vi, u)| vi * u).collect();
        let mut v7: Vec<Scalar> = vec![Scalar::zero(); 2*n+3];
        v7.extend_from_slice(&u_v_s);
        v7.extend_from_slice(&vec![Scalar::zero(); sn]);

        // theta
        let theta: Vec<Scalar> = v0.iter().zip(v1.iter()).map(|(v0i, v1i)| v0i + z * v1i).collect();
       
        // theta_inv
        let theta_inv: Vec<Scalar> = theta.iter().map(|thetai| thetai.invert()).collect();

        // zeta
        let z2 = util::scalar_exp_vartime(&z, 2);
        let z3 = z*z2;
        let z4 = util::scalar_exp_vartime(&z, 4);
        let z5 = z*z4;
        let z6 = z2*z4;
        let zeta: Vec<Scalar> = (0..t).map(|i| z2*v2[i] + z3*v3[i] + z4*v4[i] + z5*v5[i] + z6*v6[i]).collect();

        // nu
        let nu: Vec<Scalar> = (0..t).map(|i| z2*v7[i] + z6*v6[i]).collect();

        // alpha
        // Note that alpha = theta_inv * nu
        // No (-) sign before nu since we want <cL+cR-1^{t}, v6>=0
        let alpha: Vec<Scalar> = (0..t).map(|i| theta_inv[i] * nu[i]).collect();

        // beta
        let beta: Vec<Scalar> = (0..t).map(|i| theta_inv[i] * zeta[i]).collect();

        // delta
        let ones_ys: Scalar = y_s.iter().sum();
        let ones_ys_1: Scalar = ones_ys + util::scalar_exp_vartime(&y, s as u64);
        let onet_v6: Scalar = y_sn.iter().sum();
        let delta_cons: Scalar = z*ones_ys + z5*ones_ys_1 + z6*onet_v6;
        let delta: Scalar = (0..t).map(|i| alpha[i] * zeta[i]).fold(delta_cons, |acc, x| acc + x);
        
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
pub struct MProvePlus {
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

impl MProvePlus {
    pub fn prove(
        // crs
        G: &RistrettoPoint,            // base for amount in a pedersen commitment
        H: &RistrettoPoint,            // base for blinding factor in a pedersen commitment
        Gt: &RistrettoPoint,           // instead of G1
        H_prime: &RistrettoPoint,      // base used in commitment to inner product
        p_vec: &[RistrettoPoint], 
        g_prime_vec: &[RistrettoPoint],
        h_vec: &[RistrettoPoint],
        g_vec_append: &[RistrettoPoint],
        h_vec_append: &[RistrettoPoint],
        // stmt
        C_vec: &[RistrettoPoint],      // vector of commitments
        P_vec: &[RistrettoPoint],      // addresses in the ring (public keys)
        H_vec: &[RistrettoPoint],      // hash of addresses
        // wit
        E_vec: &[Scalar],              // secret indices
        x_vec: &[Scalar],              // secret keys
        gamma: &Scalar,                // blinding factor of total reserves
    ) -> MProvePlus {

        // Ring size
        let n = P_vec.len();

        // number of addresses owned by the exchange
        let s = x_vec.len();
        
        // size of honestly encoded witness vector
        let t = s*n + 2*n + s + 3;

        // transcript used for generating challenges
        let mut transcript: Vec<u8> = Vec::new();

        // other prelims
        let p_len = 2*n + s + 3;
        let mut rng = rand::thread_rng();

        // generate u,v
        transcript.extend_from_slice(G.compress().as_bytes());
        let u = Scalar::hash_from_bytes::<Sha512>(&transcript);
        transcript.extend_from_slice(H.compress().as_bytes());
        let v = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // define e_prime, E_mat and e_hat
        let e_prime = E_vec;
        let mut E_mat: Vec<Scalar> = vec![Scalar::zero(); s*n];
        
        let mut I_vec_temp: Vec<RistrettoPoint> = vec![*G; s];
        let mut C_vec_temp: Vec<RistrettoPoint> = vec![*G; s];

        let v_s: Vec<Scalar> = util::exp_iter(v).take(s).collect();
        
        let mut index: usize = 0;
        let e_hat: Vec<Scalar> = (0..n)
            .map(|i| {
                if E_vec[i] == Scalar::one() {
                    // update 'index' numbered row of E_mat
                    E_mat[n*index + i] = Scalar::one();
                    
                    // generate key-images I = (x_i)H_i
                    I_vec_temp[index] = H_vec[i] * x_vec[index];
                    // I_vec_temp.push(H_vec[i] * x_vec[index]);

                    // update secret commitment vector
                    // C_vec_temp.push(C_vec[i]);
                    C_vec_temp[index] = C_vec[i];

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
        let sn = E_mat.len();
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

        // generate key-images from I_vec_temp
        let I_vec = I_vec_temp;

        // generate commitment to total assets
        let gamma_Gt = Gt * gamma;
        let one_vec: Vec<Scalar> = (0..s).map(|_| Scalar::one()).collect();
        let C_res1 = RistrettoPoint::vartime_multiscalar_mul(one_vec.iter(), C_vec_temp.iter());
        let C_res = gamma_Gt + C_res1;

        // define compressed stmt Y_hat_vec
        let P_vec_temp: Vec<RistrettoPoint> = (0..n).map(|i| P_vec[i] * u).collect();

        let u_square = u*u;
        let minus_u_sq = -u_square;
        let H_vec_temp: Vec<RistrettoPoint> = (0..n).map(|i| H_vec[i] * u_square).collect();
        
        let Y_hat_vec: Vec<RistrettoPoint> = (0..n).map(|i| H_vec_temp[i] + P_vec_temp[i]).collect();

        // define compressed I_hat_vec
        let I_hat_vec: Vec<RistrettoPoint> = (0..s).map(|i| I_vec[i] * (minus_u_sq * v_s[i])).collect();

        // define x_inv_vec
        let x_inv_vec: Vec<Scalar> = (0..s).map(|i| x_vec[i].invert()).collect();

        // define compressed secrets
        let minus_u = -u;
        let xi = (0..s).map(|i| x_vec[i] * minus_u * v_s[i]).fold(Scalar::zero(), |acc, x| acc + x);

        // secret vectors
        // c L := ( ξ || − 1 || γ || ê || e' || x^◦−1 || vec(E))
        let mut c_L: Vec<Scalar> = vec![xi.clone()];
        c_L.push(-Scalar::one());
        c_L.push(*gamma);
        c_L.extend_from_slice(&e_hat.clone());
        c_L.extend_from_slice(&e_prime.clone());
        c_L.extend_from_slice(&x_inv_vec.clone());
        c_L.extend_from_slice(&E_mat.clone());
        
        // c R := (0^2n+3 || x || vec(E) − 1^sn )
        let mut c_R: Vec<Scalar> = vec![Scalar::zero(); 2*n+3];
        c_R.extend_from_slice(&x_vec);
        c_R.extend_from_slice(&E_mat_comp);
    
        // defining g_0
        let mut g_vec_0: Vec<RistrettoPoint> = Vec::new();
        g_vec_0.extend_from_slice(&p_vec);
        g_vec_0.extend_from_slice(&g_prime_vec);

        // P -> V: A
        // A = (H' * r_a) + <g_0 * c_L> + <h * c_R>
        let r_A = Scalar::random(&mut rng);
        
        let A = RistrettoPoint::vartime_multiscalar_mul(
            c_L.iter().chain(c_R.iter()).chain(iter::once(&r_A)),
            g_vec_0.iter().chain(h_vec.iter()).chain(iter::once(H_prime)),
        );

        // challenge w
        transcript.extend_from_slice(A.compress().as_bytes());
        let w = Scalar::hash_from_bytes::<Sha512>(&transcript);

        // defining g_w
        let mut g_vec_w: Vec<RistrettoPoint> = vec![*G, C_res, *Gt];
        g_vec_w.extend_from_slice(&Y_hat_vec);
        g_vec_w.extend_from_slice(&C_vec);
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
        let s_R = (0..t).map(|_| Scalar::random(&mut rng)).collect::<Vec<Scalar>>();
        let s_L = (0..t).map(|_| Scalar::random(&mut rng)).collect::<Vec<Scalar>>();

        let H_prime_r_S = H_prime * r_S;
        let sL_gw = (0..t).map(|i| g_vec_w[i] * s_L[i]).fold(H_prime_r_S, |acc, x| acc + x);
        let S = (0..t).map(|i| h_vec[i] * s_R[i]).fold(sL_gw, |acc, x| acc + x);

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
        
        // calculate t2, t1, t0
        let t2 = (0..t).map(|i| (s_L[i] * s_R[i]) * theta[i]).fold(Scalar::zero(), |acc, x| acc + x);
        
        let t1 = (0..t)
            .map(|i| {
                (c_L[i] + alpha[i])*theta[i]*s_R[i] + (c_R[i]*theta[i] + zeta[i])*s_L[i]
            })
            .fold(Scalar::zero(), |acc, x| acc + x); 

        // P -> V: T_1, T_2
        let tau1 = Scalar::random(&mut rng);
        let tau2 = Scalar::random(&mut rng);
        let T1 = G*t1 + H*tau1;
        let T2 = G*t2 + H*tau2;

        // generate challenge x
        transcript.extend_from_slice(T1.compress().as_bytes());
        transcript.extend_from_slice(T2.compress().as_bytes());
        let x = Scalar::hash_from_bytes::<Sha512>(&transcript);
        
        // P -> V: tau_x, r, t_hat
        // compute tau_x, r, Lp, Rp, t_hat
        let taux_1 = x * tau1;
        let taux_2 = x * x * tau2;
        let tau_x = taux_1 + taux_2;

        let r = (r_S*x) + r_A;
        
        let Lp: Vec<Scalar> = (0..t).map(|i| s_L[i]*x + (c_L[i] + alpha[i])).collect();

        let Rp: Vec<Scalar> = (0..t).map(|i| theta[i]*(c_R[i] + s_R[i]*x) + zeta[i]).collect();

        let t_hat = Lp.iter().zip(Rp.iter()).fold(Scalar::zero(), |acc, x| {
            acc + (x.0*x.1)
        });

        // Running inner product argument
        let Q = RistrettoPoint::hash_from_bytes::<Sha512>(b"test point");

        // Run ipp for non-power of two secret vectors
        let N = t.next_power_of_two();
        let res = N-t;
        let zero_append_vec = vec![Scalar::zero();res];

        // Append 0s to secret vectors
        let mut a = Lp.clone();
        let mut b = Rp.clone();
        a.extend_from_slice(&zero_append_vec);
        b.extend_from_slice(&zero_append_vec);

        // Multipliers of base vectors
        let G_factors: Vec<Scalar> = iter::repeat(Scalar::one()).take(N).collect();
        let mut H_factors: Vec<Scalar> = theta_inv.clone();
        H_factors.extend_from_slice(&vec![Scalar::one(); res]);

        // Append random group generators to g_vec_w and hi_tag        
        let mut g_vec_long = g_vec_w.clone();
        let mut h_vec_long: Vec<RistrettoPoint> = Vec::with_capacity(N);
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

        return MProvePlus {
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
        Gt: &RistrettoPoint, // instead of G1
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
        let t = s*n + 2*n + s + 3;
        let p_len = 2*n+3+s;
        let N = t.next_power_of_two();
        let res = N-t;

        // transcript initialization
        let mut transcript: Vec<u8> = Vec::new();

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
        
        // we are using already computed constraint vectors by prover
        let theta_inv = self.constraint_vec.theta_inv.clone();
        let alpha = self.constraint_vec.alpha.clone();
        let delta = self.constraint_vec.delta.clone();
        let beta = self.constraint_vec.beta.clone();     

        // verification equation #2
        // lhs
        let Gt_hat = G * self.t_hat;
        let Htau_x = H * self.tau_x;
        let left_side = Gt_hat + Htau_x;

        // rhs
        let Gdelta = G * delta;
        let Tx = self.T1 * x;
        let Tx_sq = self.T2 * (x*x);
        let right_side = Gdelta + Tx + Tx_sq;

        // towards verification eqn #3
        // define compressed stmt Y_hat_vec
        let P_vec_temp: Vec<RistrettoPoint> = (0..n).map(|i| P_vec[i] * u).collect();
        
        let u_square = u*u;
        let minus_u_sq = -u_square;
        let H_vec_temp: Vec<RistrettoPoint> = (0..n).map(|i| H_vec[i] * u_square).collect();
        
        let Y_hat_vec: Vec<RistrettoPoint> = (0..n).map(|i| H_vec_temp[i] + P_vec_temp[i]).collect();

        // define compressed I_hat_vec
        let v_s: Vec<Scalar> = util::exp_iter(v).take(s).collect();
        let I_hat_vec: Vec<RistrettoPoint> = (0..s).map(|i| self.I_vec[i] * minus_u_sq * v_s[i]).collect();  
        
        // defining g_w
        let mut g_vec_w: Vec<RistrettoPoint> = vec![*G, self.C_res, *Gt];
        g_vec_w.extend_from_slice(&Y_hat_vec);
        g_vec_w.extend_from_slice(&C_vec);
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
            alpha.iter().chain(beta.iter()).chain(iter::once(&self.t_hat)).chain(iter::once(&x)).chain(iter::once(&minus_r)).chain(iter::once(&Scalar::one())),
            g_vec_w.iter().chain(h_vec.iter()).chain(iter::once(&Q)).chain(iter::once(&self.S)).chain(iter::once(H_prime)).chain(iter::once(&self.A)),
        );
        
        // Append random group generators (same as prover's) to g_vec_w and hi_tag
        let mut g_vec_long = g_vec_w.clone();
        let mut h_vec_long: Vec<RistrettoPoint> = Vec::with_capacity(N);
        g_vec_long.extend_from_slice(&g_vec_append);
        h_vec_long.extend_from_slice(&h_vec);
        h_vec_long.extend_from_slice(&h_vec_append);

        // Multipliers of base vectors
        let G_factors: Vec<Scalar> = iter::repeat(Scalar::one()).take(N).collect();
        let mut H_factors: Vec<Scalar> = theta_inv.clone();
        H_factors.extend_from_slice(&vec![Scalar::one(); res]);

        let mut verifier = Transcript::new(b"innerproduct");
        let verify = self.inner_product_proof.verify(N, &mut verifier, &G_factors,
            &H_factors, &P, &Q, &g_vec_long, &h_vec_long);

        // check all three conditions are true
        if verify.is_ok() && left_side==right_side{
            Ok(())
        } else {
            Err(MProvePlusError)
        }
    }

    pub fn gen_params(n: usize, s: usize) -> (
        RistrettoPoint,
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
        Scalar,
    ) {
        
        let sn = s * n;
        let t = sn + 2*n + s + 3;
        let one = Scalar::one();
        let mut rng = rand::thread_rng();

        // generate random amounts in range {0,..,2^{32}-1}
        let a_vec: Vec<Scalar> = (0..s).map(|_| Scalar::from(rng.gen::<u32>())).collect();
        
        // generate blinding factors
        let r_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();

        // generate secret keys
        let x_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();
        
        // G, H, Gt - curve points for generating outputs and key-images
        let G = constants::RISTRETTO_BASEPOINT_POINT;
        let H = RistrettoPoint::hash_from_bytes::<Sha512>(G.compress().as_bytes());
        let Gt = RistrettoPoint::hash_from_bytes::<Sha512>(b"block_height_100");
        let H_prime = RistrettoPoint::hash_from_bytes::<Sha512>(b"h_prime");

        let pgens = PedersenGens{B: H, B_blinding: G};

        // generate p_vec, g_prime_vec, h_vec
        let p_len = 2*n+s+3;
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
        let idx = (0..s).map(|_| {
            
            let dist1 = Uniform::from(start_idx..end_idx);
            start_idx = setsize + start_idx;
            end_idx =  cmp::min(n-1, end_idx + setsize);

            dist1.sample(&mut rng)
        })
        .collect::<Vec<usize>>();

        // generate commitment to total reserves
        let gamma = Scalar::random(&mut rng);
        let Gt_gamma = Gt * gamma;
        let mut C_res = Gt_gamma;

        let mut index = 0;
        let E_vec = (0..n)
            .map(|i| {
                if index < idx.len() {
                    if i == idx[index] {
                        // generate commitments using a_vec, r_vec
                        C_vec_mut[i as usize] = pgens.commit(a_vec[index as usize], r_vec[index as usize]);
                        P_vec[i as usize] = G * x_vec[index];
                        C_res = C_res + C_vec_mut[i];
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

        (G, H, Gt, H_prime, p_vec, g_prime_vec, h_vec, g_vec_append, h_vec_append, C_vec_mut, P_vec, H_vec, E_vec, x_vec, gamma)
    }
}


#[cfg(test)]
mod tests {
    
    use super::*;
    use core::cmp;
    use rand::distributions::{Distribution, Uniform};
    use rand::Rng;
    use curve25519_dalek::constants;
    use crate::generators::{PedersenGens};
    use libc_print::{libc_println};
    use time::PreciseTime;

    #[test]
    pub fn test_mprove_plus_aok_1_in_8(){
        
        let n=8;
        let s=1;
        let sn = s * n;
        let t = sn + 2*n + s + 3;
        let one = Scalar::one();
        let mut rng = rand::thread_rng();

        // generate random amounts in range {0,..,2^{32}-1}
        let a_vec: Vec<Scalar> = (0..s).map(|_| Scalar::from(rng.gen::<u32>())).collect();
        
        // generate blinding factors
        let r_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();

        // generate secret keys
        let x_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();
        
        // G, H, Gt - curve points for generating outputs and key-images
        let G = constants::RISTRETTO_BASEPOINT_POINT;
        let H = RistrettoPoint::hash_from_bytes::<Sha512>(G.compress().as_bytes());
        let Gt = RistrettoPoint::hash_from_bytes::<Sha512>(b"block_height_100");
        let H_prime = RistrettoPoint::hash_from_bytes::<Sha512>(b"h_prime");

        let pgens = PedersenGens{B: H, B_blinding: G};

        // generate p_vec, g_prime_vec, h_vec
        let p_len = 2*n+s+3;
        let g_prime_len = t-p_len;
        let h_len = t;
        let N = (t as u64).next_power_of_two();
        let res = N-t;

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
        let Gt_gamma = Gt * gamma;
        let mut C_res = Gt_gamma;
        
        let idx = vec![6];
        let mut index = 0;
        let E_vec = (0..n)
            .map(|i| {
                if index < idx.len() {
                    if i == idx[index] {
                        // generate commitments using a_vec, r_vec
                        C_vec_mut[i as usize] = pgens.commit(a_vec[index as usize], r_vec[index as usize]);
                        P_vec[i as usize] = G * x_vec[index];
                        C_res = C_res + C_vec_mut[i as usize];
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

        // let mut transcript = Transcript::new(b"mprovetest");
        let mprove_test = MProvePlus::prove(&G, &H, &Gt, &H_prime, &p_vec, &g_prime_vec, &h_vec, &g_vec_append, &h_vec_append, &C_vec_mut, &P_vec, &H_vec, &E_vec, &x_vec, &gamma);
        
        // let mut transcript = Transcript::new(b"mprovetest");
        let result = mprove_test.verify(&G, &H, &Gt, &H_prime, &p_vec, &g_prime_vec, &h_vec, &g_vec_append, &h_vec_append, &C_vec_mut, &P_vec, &H_vec);
        assert!(result.is_ok());
    }

    pub fn test_mprove_plus_helper(n: usize, s: usize){
        
        let sn = s * n;
        let t = sn + 2*n + s + 3;
        let one = Scalar::one();
        let mut rng = rand::thread_rng();

        // generate random amounts in range {0,..,2^{32}-1}
        let a_vec: Vec<Scalar> = (0..s).map(|_| Scalar::from(rng.gen::<u32>())).collect();
        
        // generate blinding factors
        let r_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();

        // generate secret keys
        let x_vec: Vec<Scalar> = (0..s).map(|_| Scalar::random(&mut rng)).collect();
        
        // G, H, Gt - curve points for generating outputs and key-images
        let G = constants::RISTRETTO_BASEPOINT_POINT;
        let H = RistrettoPoint::hash_from_bytes::<Sha512>(G.compress().as_bytes());
        let Gt = RistrettoPoint::hash_from_bytes::<Sha512>(b"block_height_100");
        let H_prime = RistrettoPoint::hash_from_bytes::<Sha512>(b"h_prime");

        let pgens = PedersenGens{B: H, B_blinding: G};

        // generate p_vec, g_prime_vec, h_vec
        let p_len = 2*n+s+3;
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
        let idx = (0..s).map(|_| {
            
            let dist1 = Uniform::from(start_idx..end_idx);
            start_idx = setsize + start_idx;
            end_idx =  cmp::min(n-1, end_idx + setsize);

            dist1.sample(&mut rng)
        })
        .collect::<Vec<usize>>();

        // generate commitment to total reserves
        let gamma = Scalar::random(&mut rng);
        let Gt_gamma = Gt * gamma;
        let mut C_res = Gt_gamma;

        let mut index = 0;
        let E_vec = (0..n)
            .map(|i| {
                if index < idx.len() {
                    if i == idx[index] {
                        // generate commitments using a_vec, r_vec
                        C_vec_mut[i as usize] = pgens.commit(a_vec[index as usize], r_vec[index as usize]);
                        P_vec[i as usize] = G * x_vec[index];
                        C_res = C_res + C_vec_mut[i];
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

            libc_println!("For n={}, s={},", n, s);
            let start = PreciseTime::now();
            // let mut transcript = Transcript::new(b"mprovetest");
            let mprove_test = MProvePlus::prove(&G, &H, &Gt, &H_prime, &p_vec, &g_prime_vec, &h_vec, &g_vec_append, &h_vec_append, &C_vec_mut, &P_vec, &H_vec, &E_vec, &x_vec, &gamma);
            let end = PreciseTime::now();
            libc_println!("Generation time: {:?}", start.to(end));
            
            let start = PreciseTime::now();
            // let mut transcript = Transcript::new(b"mprovetest");
            let result = mprove_test.verify(&G, &H, &Gt, &H_prime, &p_vec, &g_prime_vec, &h_vec, &g_vec_append, &h_vec_append, &C_vec_mut, &P_vec, &H_vec);
            let end = PreciseTime::now();
            libc_println!("Verification time: {:?}", start.to(end));

            assert!(result.is_ok());
    }

    #[test]
    fn mPlus_2_in_100(){
        test_mprove_plus_helper(100, 2);
    }
    
    #[test]
    fn mPlus_5_in_100(){
        test_mprove_plus_helper(100, 5);
    }

    #[test]
    fn mPlus_4_in_500(){
        test_mprove_plus_helper(500, 4);
    }

    #[test]
    fn mPlus_19_in_20(){
        test_mprove_plus_helper(20, 19);
    }

}
