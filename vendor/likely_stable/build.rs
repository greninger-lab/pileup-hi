use rustc_version::{version, version_meta, Channel,Version};

fn main() {

    if version().unwrap() >= Version::parse("1.52.0").unwrap() {
        println!("crago::rustc-cfg=const_checked_div")
    }

    if let Channel::Nightly = version_meta().unwrap().channel {
        println!("crago::rustc-cfg=nightly")
    }
}
