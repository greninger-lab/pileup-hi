build:
	# cargo build --release && cp ./target/release/indel_left_align ~/.cargo/bin/indel
	cargo build --release && cp ./target/release/pileuphi ~/.cargo/bin/pileuphi

build_test:
	cargo build && cp ./target/debug/pileuphi ~/.cargo/bin/pileuphi

test:
	cargo test -- --show-output

