# Plot 1			
cargo build --release
cargo run --release --bin omnires_ristretto_bin 150 100 -n 1
cargo run --release --bin omnires_ristretto_bin 250 100 -n 1
cargo run --release --bin omnires_ristretto_bin 300 100 -n 1
cargo run --release --bin omnires_ristretto_bin 400 100 -n 1
cargo run --release --bin omnires_ristretto_bin 600 100 -n 1
cargo run --release --bin omnires_ristretto_bin 1000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 1200 100 -n 1
cargo run --release --bin omnires_ristretto_bin 2000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 2500 100 -n 1
cargo run --release --bin omnires_ristretto_bin 4500 100 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 8000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 10000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 16000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 30000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 35000 100 -n 1
# cargo run --release --bin omnires_ristretto_bin 60000 100 -n 1
# cargo run --release --bin omnires_ristretto_bin 80000 100 -n 1

# Plot 2
cargo run --release --bin omnires_ristretto_bin 5000 50 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 50 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 90 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 200 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 500 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 800 -n 1
# cargo run --release --bin omnires_ristretto_bin 5000 1200 -n 1
# cargo run --release --bin omnires_ristretto_bin 5000 1500 -n 1
# cargo run --release --bin omnires_ristretto_bin 5000 2000 -n 1
# cargo run --release --bin omnires_ristretto_bin 5000 2500 -n 1

# Heavy memory usage simulations to be left later to run
cargo run --release --bin omnires_ristretto_bin 60000 100 -n 1
cargo run --release --bin omnires_ristretto_bin 80000 100 -n 1

cargo run --release --bin omnires_ristretto_bin 5000 1200 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 1500 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 2000 -n 1
cargo run --release --bin omnires_ristretto_bin 5000 2500 -n 1
