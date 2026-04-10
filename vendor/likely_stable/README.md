[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE-MIT)
[![LICENSE](https://img.shields.io/badge/license-apache-blue.svg)](LICENSE-APACHE)
[![Documentation](https://docs.rs/likely_stable/badge.svg)](https://docs.rs/likely_stable)
[![Crates.io Version](https://img.shields.io/crates/v/likely_stable.svg)](https://crates.io/crates/likely_stable)

This crates brings [likely](https://doc.rust-lang.org/core/intrinsics/fn.likely.html) 
and [unlikely](https://doc.rust-lang.org/core/intrinsics/fn.unlikely.html) branch prediction hints to stable rust
```rust
use likely_stable::{likely,unlikely};
use rand::random;

if likely(random::<i32>() > 10) {
    println!("likely!")
} else {
    println!("unlikely!")
}
```

It also provides `if_likely` and `if_unlikely` for branch prediction
for `if let` statements.
```rust
use likely_stable::if_likely;
use rand::random;

let v = Some(random()).filter(|v:&i32| *v > 10);

if_likely!{let Some(v) = v => {
    println!("likely!")
} else {
    println!("unlikely!")
}};
```

Moreover traits `LikelyBool`, `LikelyOption` and `LikelyResult` provides *likely*
and *unlikely* versions of the methods commonly used for types `bool`, `Option` and
`Result`
```rust
use likely_stable::LikelyOption;
use rand::random;

let v = Some(random()).filter(|v:&i32| *v > 10);

v.map_or_else_likely(
    || println!("unlikely"),
    |v| println!("likely {}",v));
```

# Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
likely_stable = "0.1"
```
