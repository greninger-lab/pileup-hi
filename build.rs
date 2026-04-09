fn main() {
    cc::Build::new()
        .file("src/overlap_wrapper.c")
        .includes(vec!["./vendor/hts-sys/"])
        .includes(vec!["./vendor/hts-sys/htslib/"])
        // pkg_config::Config::new().probe("htslib")
        // .expect("failed to find htslib installation on this system").include_paths
        // )
        .opt_level(3)
        .flag_if_supported("-march=native")
        .compile("overlap_wrapper");
}
