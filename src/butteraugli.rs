use image::{DynamicImage, GenericImageView};
use std::path::PathBuf;
use std::process;

pub fn open_image(path: PathBuf) -> DynamicImage {
    let image_open_result = image::open(path);

    let img = match image_open_result {
        Ok(file) => file,
        Err(error) => {
            eprintln!("{}", error);
            process::exit(1);
        }
    };

    // The dimensions method returns the images width and height.
    println!("dimensions {:?}", img.dimensions());

    println!("{:?}", img.color());

    // The color method returns the image's `ColorType`.
    println!("{:?}", img.color());

    img
}

fn butteraugli_fuzzy_class(score: f64) -> f64 {
    let fuzzy_width_up = 6.07887388532;
    let fuzzy_width_down = 5.50793514384;
    let m0 = 2.0;
    let scaler = 0.840253347958;

    let mut val = 0.0;

    if score < 1.0 {
        val = m0 / (1.0 + ((score - 1.0) * fuzzy_width_down).exp());
        val -= 1.0;
        val += 2.0 - scaler;
        val += scaler;
    } else {
        val = m0 / (1.0 + ((score - 1.0) * fuzzy_width_up).exp());
        val += scaler;
    }

    val
}

#[inline]
fn dot_product(u: &[f64; 3], v: &[f64; 3]) -> f64 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

fn computer_kernel(sigma: f64) -> Vec<f64> {
    let m = 2.25;
    let scaler = -1.0 / (2.0 * sigma * sigma);
    let diff = (m * sigma.abs()).max(1.0) as isize;

    let mut kernel = Vec::with_capacity((2 * diff + 1) as usize);

    for i in -diff..=diff {
        kernel[(i + diff) as usize] = (scaler * i as f64 * i as f64).exp();
    }

    kernel
}
