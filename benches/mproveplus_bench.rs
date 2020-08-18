#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;
extern crate mProvePlus_ristretto;

mod bench_mproveplus {
    use mProvePlus_ristretto::proofs::mprove_plus_ristretto::MProvePlus;
    use criterion::Criterion;

    static OWN_SET_SIZES: [usize; 5] = [5, 10, 20, 50, 80];

    fn create_mproveplus_helper(n: usize, c: &mut Criterion) {
        let label = format!("Create MProve+ proofs with n={} and s=", n);

        c.bench_function_over_inputs(
            &label,
            move |b, &&s| {

                // create crs
                let (G, H, Gt, 
                    H_prime, 
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
                    gamma) 
                    = MProvePlus::gen_params(n, s);

                
                b.iter(|| {
                    MProvePlus::prove(&G, &H, &Gt, 
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
                        &gamma);
                })
            },
            &OWN_SET_SIZES,
        );
    }

    fn verify_mproveplus_helper(n: usize, c: &mut Criterion) {
        let label = format!("Verify MProve+ proofs with n={} and s=", n);

        c.bench_function_over_inputs(
            &label,
            move |b, &&s| {

                // create crs
                let (G, H, Gt, 
                    H_prime, 
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
                    gamma) 
                    = MProvePlus::gen_params(n, s);

                let mprove_plus_proof = MProvePlus::prove(&G, &H, &Gt, 
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
                    &gamma);
                
                b.iter(|| {
                    let result = mprove_plus_proof.verify(&G, &H, &Gt, 
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
                })
            },
            &OWN_SET_SIZES,
        );
    }

    fn fast_verify_mproveplus_helper(n: usize, c: &mut Criterion) {
        let label = format!("Fast Verify MProve+ proofs with n={} and s=", n);

        c.bench_function_over_inputs(
            &label,
            move |b, &&s| {

                // create crs
                let (G, H, Gt, 
                    H_prime, 
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
                    gamma) 
                    = MProvePlus::gen_params(n, s);

                let mprove_plus_proof = MProvePlus::prove(&G, &H, &Gt, 
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
                    &gamma);
                
                b.iter(|| {
                    let result = mprove_plus_proof.fast_verify(&G, &H, &Gt, 
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
                })
            },
            &OWN_SET_SIZES,
        );
    }

    pub fn create_mproveplus_100(c: &mut Criterion) {
        create_mproveplus_helper(100, c);
    }

    pub fn verify_mproveplus_100(c: &mut Criterion) {
        verify_mproveplus_helper(100, c);
    }

    pub fn fast_verify_mproveplus_100(c: &mut Criterion) {
        fast_verify_mproveplus_helper(100, c);
    }

    criterion_group! {
    name = mproveplus;
    config = Criterion::default().sample_size(10);
    targets =
        create_mproveplus_100,
        verify_mproveplus_100,
        fast_verify_mproveplus_100,
    }

}

//fn main() {}
criterion_main!(
    bench_mproveplus::mproveplus,
);

