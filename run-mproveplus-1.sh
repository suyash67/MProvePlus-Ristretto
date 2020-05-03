# Plot 1		                            MProvePlus		
# -----------------------------------------------------------------------------------------------
# n	        s	t = sn + 2n + s +3	 N = t.next_pow_2()	Est. gen time (sec)	Est. ver time (sec)
# -----------------------------------------------------------------------------------------------
# 250	    200	50703	             16	                20.300	             2.990
# 500	    200	101203	             17	                40.519	             5.968
# 1000	    200	202203	             18	                80.956	             11.924
# 2500	    200	505203	             19	                202.269	             29.792
# 5000	    200	1010203	             20	                404.456	             59.573
# 10000	    200	2020203	             21	                808.830	             119.133
# 25000	    200	5050203	             23	                2021.954	         297.815
# 50000	    200	10100203	         24	                4043.826	         595.618
# 100000	200	20200203	         25	                8087.571	         1191.224
# 250000	200	50500203	         26	                20218.806	         2978.041
# -----------------------------------------------------------------------------------------------
# Total (hrs) ~=	11.4			                    9.9804	             1.4700

# Note: timings estimated based on Intel® Core™ i7-5500U CPU @ 2.40GHz (on a single core)

cargo build --release
cargo run --release --bin mprove_plus_ristretto_bin 250 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 500 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 1000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 2500 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 5000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 10000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 25000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 50000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 100000 200 -n 1
cargo run --release --bin mprove_plus_ristretto_bin 250000 200 -n 1

