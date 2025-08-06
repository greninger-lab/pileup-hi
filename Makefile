build:
	cargo build --release && cp ./target/release/indel_left_align ~/.cargo/bin/indel

test:
	cargo test -- --show-output

