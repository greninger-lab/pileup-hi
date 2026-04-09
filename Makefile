build_test:
	cargo build && cp ./target/debug/pileuphi ~/.cargo/bin/pileuphi

build:
	cargo build --release && cp ./target/release/pileuphi ~/.cargo/bin/pileuphi

test:
	cargo test -- --show-output

