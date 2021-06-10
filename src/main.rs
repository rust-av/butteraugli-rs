use std::path::PathBuf;
use structopt::StructOpt;

use std::env;
use std::ffi::CString;
use std::os::raw::c_char;
use std::os::unix::ffi::OsStrExt;
use std::ptr;

#[cxx::bridge]
mod ffi {
    extern "C++" {
        include!("butteraugli-rs/butteraugli-cc/butteraugli.h");

        unsafe fn Run(argc: i32, argv: *mut *mut c_char) -> i32;
    }
}

#[derive(Debug, StructOpt)]
#[structopt(
    name = "butteraugli-rs",
    about = "Rust port of the Google butteraugli metric"
)]
struct CommandLineArgs {
    /// Usage: {*.(png|jpg|jpeg)} {*.(png|jpg|jpeg)}
    #[structopt(parse(from_os_str), required = true, number_of_values = 2)]
    input: Vec<PathBuf>,

    /// In PNM image format {*.(pnm)}
    #[structopt(parse(from_os_str), short, long)]
    output: Option<PathBuf>,
}

fn main() {
    // let args = CommandLineArgs::from_args();

    let args: Vec<CString> = env::args_os()
        .map(|os_str| {
            let bytes = os_str.as_bytes();
            CString::new(bytes).unwrap_or_else(|nul_error| {
                let nul_position = nul_error.nul_position();
                let mut bytes = nul_error.into_vec();
                bytes.truncate(nul_position);
                CString::new(bytes).unwrap()
            })
        })
        .collect();

    let argc = args.len();
    let mut argv: Vec<*mut c_char> = Vec::with_capacity(argc + 1);
    for arg in &args {
        argv.push(arg.as_ptr() as *mut c_char);
    }
    argv.push(ptr::null_mut()); // Nul terminator.

    unsafe {
        ffi::Run(argc as i32, argv.as_mut_ptr());
    }
}
