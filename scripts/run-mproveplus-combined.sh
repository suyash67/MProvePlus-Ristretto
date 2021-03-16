# Plot 1		MProvePlus			
# -----------------------------------------------------------------------------------------------
# n	    s	t = sn+2n+s+3	N = t.next_pow_2()	Est. gen time (sec)	Est. ver time (sec)
# -----------------------------------------------------------------------------------------------
# 150	100	15403	        16384	            5.429	            0.846
# 200	100	20503	        32768	            10.858	            1.692
# 300	100	30703	        32768	            10.858	            1.692
# 400	100	40903	        65536	            21.716	            3.384
# 600	100	61303	        65536	            21.716	            3.384
# 1000	100	102103	        131072	            43.432	            6.768
# 1200	100	122503	        131072	            43.432	            6.768
# 2000	100	204103	        262144	            86.864	            13.536
# 2500	100	255103	        262144	            86.864	            13.536
# 4500	100	459103	        524288	            173.728	            27.072
# 5000	100	510103	        524288	            173.728	            27.072
# 8000	100	816103	        1048576	            347.456	            54.144
# 10000	100	1020103	        1048576	            347.456	            54.144
# 16000	100	1632103	        2097152	            694.912	            108.288
# 30000	100	3060103	        4194304	            1389.824	        216.576
# 35000	100	3570103	        4194304	            1389.824	        216.576
# 60000	100	6120103	        8388608	            2779.648	        433.152
# 80000	100	8160103	        8388608	            2779.648	        433.152
# -----------------------------------------------------------------------------------------------
# Total (hrs) =	3.4452			                2.8894          	0.4503

# Note: timings ESTIMATED based on Intel® Core™ i7-5500U CPU @ 2.40GHz (on a single core)

cargo build --release
cargo run --release --bin mprove_plus_ristretto_bin 150 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 250 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 300 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 400 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 600 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 1000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 1200 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 2000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 2500 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 4500 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 8000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 16000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 30000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 35000 100 -n 1
# cargo run --release --bin mprove_plus_ristretto_bin 60000 100 -n 1
# cargo run --release --bin mprove_plus_ristretto_bin 80000 100 -n 1

# Plot 2		MProvePlus	
# -----------------------------------------------------------------------------------------------
# n	    s	      t = sn+2n+s+3	 N = t.next_pow_2()	Est. gen time (sec)	Est. ver time (sec)
# -----------------------------------------------------------------------------------------------
# 5000	45	      235048	     262144	            100.24	            16.456
# 5000	50	      260053	     262144	            100.24	            16.456
# 5000	90	      460093	     524288	            200.48	            32.912
# 5000	100	      510103	     524288	            200.48	            32.912
# 5000	200	      1010203	     1048576            400.96              65.824
# 5000	500	      2510503	     4194304            1603.84             263.296
# 5000	800	      4010803	     4194304            1603.84             263.296
# 5000	1200	  6011203	     8388608            3207.68             526.592
# 5000	1500	  7511503	     8388608            3207.68             526.592
# 5000	2000	  10012003	     16777216           6415.36             1053.184
# 5000	2500	  12512503	     16777216           6415.36             1053.184
# -----------------------------------------------------------------------------------------------
# Total (hrs) =	6.5823			6.5156	0.0503

cargo run --release --bin mprove_plus_ristretto_bin 5000 50 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 50 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 90 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 500 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 800 -n 1
# cargo run --release --bin mprove_plus_ristretto_bin 5000 1200 -n 1
# cargo run --release --bin mprove_plus_ristretto_bin 5000 1500 -n 1
# cargo run --release --bin mprove_plus_ristretto_bin 5000 2000 -n 1
# cargo run --release --bin mprove_plus_ristretto_bin 5000 2500 -n 1

# Heavy memory usage simulations to be left later to run
cargo run --release --bin mprove_plus_ristretto_bin 60000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 80000 100 -n 1

cargo run --release --bin mprove_plus_ristretto_bin 5000 1200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 1500 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 2000 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 2500 -n 1
