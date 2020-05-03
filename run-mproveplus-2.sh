# Plot 2		MProvePlus	
# -----------------------------------------------------------------------------------------------
# n	    s	    t = sn + 2n + s +3	N = t.next_pow_2()	Est. gen time (sec)	Est. ver time (sec)
# -----------------------------------------------------------------------------------------------
# 10000	50	    520053	            19	                196.792	            32.128
# 10000	100	    1020103	            20	                386.0147131	        63.02024829
# 10000	200	    2020203	            21	                764.4601392	        124.8047449
# 10000	500	    5020503	            23	                1899.796418	        310.1582346
# 10000	1000	10021003	        24	                3792.023548	        619.0807175
# 10000	2000	20022003	        25	                7576.47781	        1236.925683
# 10000	5000	50025003	        26	                18929.84059	        3090.460581
# -----------------------------------------------------------------------------------------------
# Total (hrs) ~= 10.8			                        9.3182	            1.5213

# Note: timings estimated based on Intel® Core™ i7-5500U CPU @ 2.40GHz (on a single core)

cargo build --release
cargo run --release --bin mprove_plus_ristretto_bin 10000 50 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 100 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 500 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 1000 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 2000 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 5000 -n 1