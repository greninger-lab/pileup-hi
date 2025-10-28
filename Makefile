build:
	# cargo build --release && cp ./target/release/indel_left_align ~/.cargo/bin/indel
	cargo build --release && cp ./target/release/viggo ~/.cargo/bin/viggo

build_test:
	cargo build && cp ./target/debug/viggo ~/.cargo/bin/viggo

test:
	cargo test -- --show-output

