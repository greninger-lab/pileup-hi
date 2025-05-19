let 
  rust_overlay = import (builtins.fetchTarball "https://github.com/oxalica/rust-overlay/archive/master.tar.gz");

pkgs = import <nixpkgs> { overlays = [ rust_overlay ]; };
in
pkgs.mkShell {
  packages = [ pkgs.htslib pkgs.clang pkgs.rust-bin.nightly.latest.default ];
  LIBCLANG_PATH = "${pkgs.llvmPackages_16.libclang.lib}/lib";
}
