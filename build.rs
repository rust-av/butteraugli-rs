fn main() {
    cxx_build::bridge("src/main.rs")
        .file("butteraugli-cc/butteraugli.cc")
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-O3")
        .compile("cxx");
}
