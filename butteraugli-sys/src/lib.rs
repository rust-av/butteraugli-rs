/// ```
/// extern crate butteraugli_sys;
///
/// fn main() {
///    butteraugli_sys::run();
///
///    println!(
///        "{}",
///        butteraugli_sys::butteraugli(
///            "test-images/1a.png".to_owned(),
///            "test-images/1b.png".to_owned()
///        )
///    );
///}
/// ```
use std::env;
use std::ffi::CString;
use std::os::raw::c_char;
use std::os::unix::ffi::OsStrExt;
use std::ptr;

#[cxx::bridge]
mod ffi {
    extern "C++" {
        include!("butteraugli-sys/butteraugli-cc/butteraugli.h");

        unsafe fn Run(argc: i32, argv: *mut *mut c_char) -> i32;

        unsafe fn diff_value(filename1: *mut c_char, filename2: *mut c_char) -> f64;
    }
}

pub fn run() {
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

pub fn butteraugli(filename1: String, filename2: String) -> f64 {
    unsafe {
        ffi::diff_value(
            CString::new(filename1).expect("").as_ptr() as *mut c_char,
            CString::new(filename2).expect("").as_ptr() as *mut c_char,
        )
    }
}
