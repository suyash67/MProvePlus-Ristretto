#![allow(non_snake_case)]

extern crate structopt;
extern crate mProvePlus_ristretto;
extern crate libc_print;
extern crate time;

use structopt::StructOpt;
use time::{PreciseTime};
use mProvePlus_ristretto::proofs::mprove_plus_ristretto::MProvePlus;
use libc_print::{libc_println};

#[derive(Debug, StructOpt)]
#[structopt(name = "mprove-plus", about = "MProvePlus proof generation simulator using Ristretto.")]
struct Opt {
  //#[structopt(short = "a", long = "anonsize")]
  anon_list_size: usize,
  //#[structopt(short = "o", long = "ownsize")]
  own_list_size: usize,
  #[structopt(short = "n", long = "numiter", default_value = "1")]
  num_iter: u32,
}

fn main() {
    // 
    // cargo run --release --bin mprove_plus_ristretto_bin 1000 100 -n 10
    //
    let opt = Opt::from_args();

    let num_iter = opt.num_iter;
    let mut gen_proof_start;
    let mut gen_proof_end;
    let mut ver_proof_start;
    let mut ver_proof_end;
    let mut total_gen_proof_duration: i64 = 0;
    let mut total_ver_proof_duration: i64 = 0;

    let (G, H, Gt, H_prime, p_vec, g_prime_vec, h_vec, g_vec_append, h_vec_append, C_vec_mut, P_vec, H_vec, E_vec, x_vec, gamma) = MProvePlus::gen_params(opt.anon_list_size, opt.own_list_size);

    let sim_start = PreciseTime::now();

    for _i in 0..num_iter {        

        gen_proof_start = PreciseTime::now();
        let mprove_plus_proof = MProvePlus::prove(&G, &H, &Gt, &H_prime, &p_vec, &g_prime_vec, &h_vec, &g_vec_append, &h_vec_append, &C_vec_mut, &P_vec, &H_vec, &E_vec, &x_vec, &gamma);
        gen_proof_end = PreciseTime::now();
        total_gen_proof_duration += (gen_proof_start.to(gen_proof_end)).num_milliseconds();
  
        ver_proof_start = PreciseTime::now();
        let result = mprove_plus_proof.verify(&G, &H, &Gt, &H_prime, &p_vec, &g_prime_vec, &h_vec, &g_vec_append, &h_vec_append, &C_vec_mut, &P_vec, &H_vec);
        assert!(result.is_ok());
        ver_proof_end = PreciseTime::now();
        total_ver_proof_duration += (ver_proof_start.to(ver_proof_end)).num_milliseconds();
      }
  
      let sim_end = PreciseTime::now();
      libc_println!("Total simulation time = {:?}", sim_start.to(sim_end));
  
      libc_println!("Options = {:?}", opt);
      libc_println!("Average proof generation time = {:?}",
        (total_gen_proof_duration as f64/(1000.0*num_iter as f64)));
      libc_println!("Average proof verification time = {:?}",
        (total_ver_proof_duration as f64/(1000.0*num_iter as f64)));

}
