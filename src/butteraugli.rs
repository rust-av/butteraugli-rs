use image::{DynamicImage, GenericImageView};
use std::ops::{Add, Deref, Mul, Neg};
use std::path::PathBuf;
use std::{mem, process};

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

// Purpose of kInternalGoodQualityThreshold:
// Normalize 'ok' image degradation to 1.0 across different versions of
// butteraugli.
static kInternalGoodQualityThreshold: f64 = 20.35;
static kGlobalScale: f64 = 1.0 / kInternalGoodQualityThreshold;

#[inline]
pub fn dot_product(u: &[f64; 3], v: &[f64; 3]) -> f64 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

pub fn computer_kernel(sigma: f64) -> Vec<f64> {
    let m = 2.25;
    let scaler = -1.0 / (2.0 * sigma * sigma);
    let diff = (m * sigma.abs()).max(1.0) as isize;

    let mut kernel = Vec::with_capacity((2 * diff + 1) as usize);

    for i in -diff..=diff {
        kernel[(i + diff) as usize] = (scaler * i as f64 * i as f64).exp();
    }

    kernel
}

#[inline]
pub fn interpolate_clamp_negative(array: &[f64], size: isize, mut ix: f64) -> f64 {
    if ix < 0.0 {
        ix = 0.0;
    }

    let baseix = ix as isize;
    let mut res = 0.0;

    if baseix >= (size - 1) {
        res = array[size as usize - 1];
    } else {
        let mix = ix - baseix as f64;
        let nextix = baseix + 1;
        res = array[baseix as usize] + mix * (array[nextix as usize] - array[baseix as usize]);
    }

    res
}

#[inline]
pub fn opsin_absorbance<K>(in0: &K, in1: &K, in2: &K, out0: &mut K, out1: &mut K, out2: &mut K)
where
    K: Mul<f64, Output = K> + Add<K, Output = K> + Add<f64, Output = K> + Copy,
    f64: Mul<K, Output = K> + Mul<f64, Output = f64>,
{
    let mixi0 = 0.254462330846;
    let mixi1 = 0.488238255095;
    let mixi2 = 0.0635278003854;
    let mixi3 = 1.01681026909;
    let mixi4 = 0.195214015766;
    let mixi5 = 0.568019861857;
    let mixi6 = 0.0860755536007;
    let mixi7 = 1.1510118369;
    let mixi8 = 0.07374607900105684;
    let mixi9 = 0.06142425304154509;
    let mixi10 = 0.24416850520714256;
    let mixi11 = 1.20481945273;

    // what is max1, max2, etc ????
    // https://github.com/google/butteraugli/blob/71b18b636b9c7d1ae0c1d3730b85b3c127eb4511/butteraugli/butteraugli.h#L501

    *out0 = mixi0 * (*in0) + mixi1 * (*in1) + mixi2 * (*in2) + mixi3;
    *out1 = mixi4 * (*in0) + mixi5 * (*in1) + mixi6 * (*in2) + mixi7;
    *out2 = mixi8 * (*in0) + mixi9 * (*in1) + mixi10 * (*in2) + mixi11;
}

pub fn gamma_min() -> f64 {
    let mut out0 = 0.0;
    let mut out1 = 0.0;
    let mut out2 = 0.0;

    opsin_absorbance(&0.0, &0.0, &0.0, &mut out0, &mut out1, &mut out2);

    out0.min(out1.min(out2))
}

pub fn gamma_max() -> f64 {
    let mut out0 = 0.0;
    let mut out1 = 0.0;
    let mut out2 = 0.0;

    opsin_absorbance(&255.0, &255.0, &255.0, &mut out0, &mut out1, &mut out2);

    out0.max(out1.max(out2))
}

pub fn simple_gamma(mut v: f64) -> f64 {
    let kGamma = 0.372322653176;
    let limit = 37.8000499603;

    let bright = v - limit;
    if bright >= 0.0 {
        let mul = 0.0950819040934;
        v -= bright * mul;
    }
    {
        let limit2 = 74.6154406429;
        let bright2 = v - limit2;
        if bright2 >= 0.0 {
            let mul = 0.01;
            v -= bright2 * mul;
        }
    }
    {
        let limit2 = 82.8505938033;
        let bright2 = v - limit2;
        if bright2 >= 0.0 {
            let mul = 0.0316722592629;
            v -= bright2 * mul;
        }
    }
    {
        let limit2 = 92.8505938033;
        let bright2 = v - limit2;
        if bright2 >= 0.0 {
            let mul = 0.221249885752;
            v -= bright2 * mul;
        }
    }
    {
        let limit2 = 102.8505938033;
        let bright2 = v - limit2;
        if bright2 >= 0.0 {
            let mul = 0.0402547853939;
            v -= bright2 * mul;
        }
    }
    {
        let limit2 = 112.8505938033;
        let bright2 = v - limit2;
        if bright2 >= 0.0 {
            let mul = 0.021471798711500003;
            v -= bright2 * mul;
        }
    }
    let offset = 0.106544447664;
    let scale = 10.7950943969;

    scale * (offset + f64::powf(v, kGamma))
}

pub struct RationalPolynomial {
    // Domain of the polynomials; they are undefined elsewhere.
    min_value: f64,
    max_value: f64,

    // Coefficients of T_n (Chebyshev polynomials of the first kind).
    // Degree 5/5 is a compromise between accuracy (0.1%) and numerical stability.
    p: [f64; 5 + 1],
    q: [f64; 5 + 1],
}

pub fn clenshaw_recursion(x: f64, coefficients: &[f64], b1: &mut f64, b2: &mut f64, index: usize) {
    if index == 0 {
        let x_b1 = x * (*b1);
        *b1 = x_b1 - (*b2) + coefficients[0];
    } else {
        let x_b1 = x * (*b1);
        let t = (x_b1 + x_b1) - (*b2) + coefficients[index];
        *b2 = *b1;
        *b1 = t;

        clenshaw_recursion(x, &coefficients, b1, b2, index - 1);
    }
}

impl RationalPolynomial {
    pub fn evaluate_polynomial(x: f64, coefficients: &mut [f64], n: usize) -> f64 {
        let mut b1 = 0.0;
        let mut b2 = 0.0;

        clenshaw_recursion(x, coefficients, &mut b1, &mut b2, n);

        b1
    }

    #[inline]
    pub fn evaluate(&mut self, x: f64) -> f64 {
        let x01 = (x - self.min_value) / (self.max_value - self.min_value);
        let xc = 2.0 * x01 - 1.0;

        let yp = RationalPolynomial::evaluate_polynomial(xc, &mut self.p, 6);
        let yq = RationalPolynomial::evaluate_polynomial(xc, &mut self.q, 6);

        if yq == 0.0 {
            0.0
        } else {
            yp / yq
        }
    }
}

#[inline]
pub fn gamma_polynomial(value: f64) -> f64 {
    let mut r = RationalPolynomial {
        min_value: 0.971783,
        max_value: 590.188894,
        p: [
            98.7821300963361,
            164.273222212631,
            92.948112871376,
            33.8165311212688,
            6.91626704983562,
            0.556380877028234,
        ],
        q: [
            1.0,
            1.64339473427892,
            0.89392405219969,
            0.298947051776379,
            0.0507146002577288,
            0.00226495093949756,
        ],
    };

    r.evaluate(value)
}

#[inline]
pub fn gamma(v: f64) -> f64 {
    gamma_polynomial(v)
}

#[inline]
pub fn remove_range_around_zero(w: f64, x: f64) -> f64 {
    if x > w {
        x - w
    } else {
        if x < w.neg() {
            x + w
        } else {
            0.0
        }
    }
}

#[inline]
pub fn amplify_range_around_zero(w: f64, x: f64) -> f64 {
    if x > w {
        x + w
    } else {
        if x < w.neg() {
            x - w
        } else {
            2.0 * x
        }
    }
}

#[inline]
pub fn xyb_low_freq_to_vals<K>(x: &K, y: &K, b_arg: &K, valx: &mut K, valy: &mut K, valb: &mut K)
where
    K: Mul<f64, Output = K> + Add<K, Output = K> + Add<f64, Output = K> + Copy,
    f64: Mul<K, Output = K> + Mul<f64, Output = f64>,
{
    let xmuli = 5.57547552483;
    let ymuli = 1.20828034498;
    let bmuli = 6.08319517575;
    let y_to_b_muli = -0.628811683685;

    /* what is this
    const V xmul(xmuli);
    const V ymul(ymuli);
    const V bmul(bmuli);
    const V y_to_b_mul(y_to_b_muli);
    */

    let b = (*b_arg) + y_to_b_muli * (*y);
    *valb = b * bmuli;
    *valx = (*x) * xmuli;
    *valy = (*y) * ymuli;
}

pub fn suppress_hf_in_brightness(hf: f64, brightness: f64, mul: f64, reg: f64) -> f64 {
    let scalar = mul * reg / (reg + brightness);
    scalar * hf
}

pub fn suppress_uf_in_brightness(hf: f64, brightness: f64, mul: f64, reg: f64) -> f64 {
    let scalar = mul * reg / (reg + brightness);
    scalar * hf
}

pub fn maximum_clamp(mut v: f64, maxval: f64) -> f64 {
    let kMul = 0.688059627878;
    if v >= maxval {
        v -= maxval;
        v *= kMul;
        v += maxval;
    } else if v < -maxval {
        v += maxval;
        v *= kMul;
        v -= maxval;
    }

    v
}

pub fn malta_unit_lf(d: &[f64], xs: isize) -> f64 {
    let xs3 = 3 * xs;
    let mut retval = 0.0;

    {
        // x grows, y constant
        let sum = d[(-4 + 4) as usize]
            + d[(-2 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(2 + 4) as usize]
            + d[(4 + 4) as usize];
        retval += sum * sum;
    }
    {
        // y grows, x constant
        let sum = d[(-xs3 - xs + 4) as usize]
            + d[(-xs - xs + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs + 4) as usize]
            + d[(xs3 + xs + 4) as usize];
        retval += sum * sum;
    }
    {
        // both grow
        let sum = d[(-xs3 - 3 + 4) as usize]
            + d[(-xs - xs - 2 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs + 2 + 4) as usize]
            + d[(xs3 + 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        // y grows, x shrinks
        let sum = d[(-xs3 + 3 + 4) as usize]
            + d[(-xs - xs + 2 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs - 2 + 4) as usize]
            + d[(xs3 - 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        // y grows -4 to 4, x shrinks 1 -> -1
        let sum = d[(-xs3 - xs + 1 + 4) as usize]
            + d[(-xs - xs + 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs - 1 + 4) as usize]
            + d[(xs3 + xs - 1 + 4) as usize];
        retval += sum * sum;
    }
    {
        //  y grows -4 to 4, x grows -1 -> 1
        let sum = d[(-xs3 - xs - 1 + 4) as usize]
            + d[(-xs - xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs + 1 + 4) as usize]
            + d[(xs3 + xs + 1 + 4) as usize];
        retval += sum * sum;
    }
    {
        // x grows -4 to 4, y grows -1 to 1
        let sum = d[(-4 - xs + 4) as usize]
            + d[(-2 - xs + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(2 + xs + 4) as usize]
            + d[(4 + xs + 4) as usize];
        retval += sum * sum;
    }
    {
        // x grows -4 to 4, y shrinks 1 to -1
        let sum = d[(-4 + xs + 4) as usize]
            + d[(-2 + xs + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(2 - xs + 4) as usize]
            + d[(4 - xs + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1__*______
        2___*_____
        3_________
        4____0____
        5_________
        6_____*___
        7______*__
        8_________ */
        let sum = d[(-xs3 - 2 + 4) as usize]
            + d[(-xs - xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs + 1 + 4) as usize]
            + d[(xs3 + 2 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1______*__
        2_____*___
        3_________
        4____0____
        5_________
        6___*_____
        7__*______
        8_________ */
        let sum = d[(-xs3 + 2 + 4) as usize]
            + d[(-xs - xs + 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs - 1 + 4) as usize]
            + d[(xs3 - 2 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2_*_______
        3__*______
        4____0____
        5______*__
        6_______*_
        7_________
        8_________ */
        let sum = d[(-xs - xs - 3 + 4) as usize]
            + d[(-xs - 2 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 2 + 4) as usize]
            + d[(xs + xs + 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2_______*_
        3______*__
        4____0____
        5__*______
        6_*_______
        7_________
        8_________ */
        let sum = d[(-xs - xs + 3 + 4) as usize]
            + d[(-xs + 2 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs - 2 + 4) as usize]
            + d[(xs + xs - 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2________*
        3______*__
        4____0____
        5__*______
        6*________
        7_________
        8_________ */

        let sum = d[(xs + xs - 4 + 4) as usize]
            + d[(xs - 2 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(-xs + 2 + 4) as usize]
            + d[(-xs - xs + 4 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2*________
        3__*______
        4____0____
        5______*__
        6________*
        7_________
        8_________ */
        let sum = d[(-xs - xs - 4 + 4) as usize]
            + d[(-xs - 2 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 2 + 4) as usize]
            + d[(xs + xs + 4 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0__*______
        1_________
        2___*_____
        3_________
        4____0____
        5_________
        6_____*___
        7_________
        8______*__ */
        let sum = d[(-xs3 - xs - 2 + 4) as usize]
            + d[(-xs - xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs + 1 + 4) as usize]
            + d[(xs3 + xs + 2 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0______*__
        1_________
        2_____*___
        3_________
        4____0____
        5_________
        6___*_____
        7_________
        8__*______ */
        let sum = d[(-xs3 - xs + 2 + 4) as usize]
            + d[(-xs - xs + 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + xs - 1 + 4) as usize]
            + d[(xs3 + xs - 2 + 4) as usize];
        retval += sum * sum;
    }

    retval
}

pub fn malta_unit(d: &[f64], xs: isize) -> f64 {
    let xs3 = 3 * xs;
    let mut retval = 0.0;

    {
        // x grows, y constant
        let sum = d[(-4 + 4) as usize]
            + d[(-3 + 4) as usize]
            + d[(-2 + 4) as usize]
            + d[(-1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(1 + 4) as usize]
            + d[(2 + 4) as usize]
            + d[(3 + 4) as usize]
            + d[(4 + 4) as usize];
        retval += sum * sum;
    }
    {
        // y grows, x constant
        let sum = d[(-xs3 - xs + 4) as usize]
            + d[(-xs3 + 4) as usize]
            + d[(-xs - xs + 4) as usize]
            + d[(-xs + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 4) as usize]
            + d[(xs + xs + 4) as usize]
            + d[(xs3 + 4) as usize]
            + d[(xs3 + xs + 4) as usize];
        retval += sum * sum;
    }
    {
        // both grow
        let sum = d[(-xs3 - 3 + 4) as usize]
            + d[(-xs - xs - 2 + 4) as usize]
            + d[(-xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 1 + 4) as usize]
            + d[(xs + xs + 2 + 4) as usize]
            + d[(xs3 + 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        // y grows, x shrinks
        let sum = d[(-xs3 + 3 + 4) as usize]
            + d[(-xs - xs + 2 + 4) as usize]
            + d[(-xs + 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs - 1 + 4) as usize]
            + d[(xs + xs - 2 + 4) as usize]
            + d[(xs3 - 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        // y grows -4 to 4, x shrinks 1 -> -1
        let sum = d[(-xs3 - xs + 1 + 4) as usize]
            + d[(-xs3 + 1 + 4) as usize]
            + d[(-xs - xs + 1 + 4) as usize]
            + d[(-xs + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 4) as usize]
            + d[(xs + xs - 1 + 4) as usize]
            + d[(xs3 - 1 + 4) as usize]
            + d[(xs3 + xs - 1 + 4) as usize];
        retval += sum * sum;
    }
    {
        //  y grows -4 to 4, x grows -1 -> 1
        let sum = d[(-xs3 - xs - 1 + 4) as usize]
            + d[(-xs3 - 1 + 4) as usize]
            + d[(-xs - xs - 1 + 4) as usize]
            + d[(-xs + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 4) as usize]
            + d[(xs + xs + 1 + 4) as usize]
            + d[(xs3 + 1 + 4) as usize]
            + d[(xs3 + xs + 1 + 4) as usize];
        retval += sum * sum;
    }
    {
        // x grows -4 to 4, y grows -1 to 1
        let sum = d[(-4 - xs + 4) as usize]
            + d[(-3 - xs + 4) as usize]
            + d[(-2 - xs + 4) as usize]
            + d[(-1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(1 + 4) as usize]
            + d[(2 + xs + 4) as usize]
            + d[(3 + xs + 4) as usize]
            + d[(4 + xs + 4) as usize];
        retval += sum * sum;
    }
    {
        // x grows -4 to 4, y shrinks 1 to -1
        let sum = d[(-4 + xs + 4) as usize]
            + d[(-3 + xs + 4) as usize]
            + d[(-2 + xs + 4) as usize]
            + d[(-1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(1 + 4) as usize]
            + d[(2 - xs + 4) as usize]
            + d[(3 - xs + 4) as usize]
            + d[(4 - xs + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1__*______
        2___*_____
        3___*_____
        4____0____
        5_____*___
        6_____*___
        7______*__
        8_________ */
        let sum = d[(-xs3 - 2 + 4) as usize]
            + d[(-xs - xs - 1 + 4) as usize]
            + d[(-xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 1 + 4) as usize]
            + d[(xs + xs + 1 + 4) as usize]
            + d[(xs3 + 2 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1______*__
        2_____*___
        3_____*___
        4____0____
        5___*_____
        6___*_____
        7__*______
        8_________ */
        let sum = d[(-xs3 + 2 + 4) as usize]
            + d[(-xs - xs + 1 + 4) as usize]
            + d[(-xs + 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs - 1 + 4) as usize]
            + d[(xs + xs - 1 + 4) as usize]
            + d[(xs3 - 2 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2_*_______
        3__**_____
        4____0____
        5_____**__
        6_______*_
        7_________
        8_________ */
        let sum = d[(-xs - xs - 3 + 4) as usize]
            + d[(-xs - 2 + 4) as usize]
            + d[(-xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 1 + 4) as usize]
            + d[(xs + 2 + 4) as usize]
            + d[(xs + xs + 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2_______*_
        3_____**__
        4____0____
        5__**_____
        6_*_______
        7_________
        8_________ */
        let sum = d[(-xs - xs + 3 + 4) as usize]
            + d[(-xs + 2 + 4) as usize]
            + d[(-xs + 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs - 1 + 4) as usize]
            + d[(xs - 2 + 4) as usize]
            + d[(xs + xs - 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2_________
        3______**_
        4____0*___
        5__**_____
        6**_______
        7_________
        8_________ */

        let sum = d[(xs + xs - 4 + 4) as usize]
            + d[(xs + xs - 3 + 4) as usize]
            + d[(xs - 2 + 4) as usize]
            + d[(xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(1 + 4) as usize]
            + d[(-xs + 2 + 4) as usize]
            + d[(-xs + 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0_________
        1_________
        2**_______
        3__**_____
        4____0*___
        5______**_
        6_________
        7_________
        8_________ */
        let sum = d[(-xs - xs - 4 + 4) as usize]
            + d[(-xs - xs - 3 + 4) as usize]
            + d[(-xs - 2 + 4) as usize]
            + d[(-xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(1 + 4) as usize]
            + d[(xs + 2 + 4) as usize]
            + d[(xs + 3 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0__*______
        1__*______
        2___*_____
        3___*_____
        4____0____
        5____*____
        6_____*___
        7_____*___
        8_________ */
        let sum = d[(-xs3 - xs - 2 + 4) as usize]
            + d[(-xs3 - 2 + 4) as usize]
            + d[(-xs - xs - 1 + 4) as usize]
            + d[(-xs - 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 4) as usize]
            + d[(xs + xs + 1 + 4) as usize]
            + d[(xs3 + 1 + 4) as usize];
        retval += sum * sum;
    }
    {
        /* 0______*__
        1______*__
        2_____*___
        3_____*___
        4____0____
        5____*____
        6___*_____
        7___*_____
        8_________ */
        let sum = d[(-xs3 - xs + 2 + 4) as usize]
            + d[(-xs3 + 2 + 4) as usize]
            + d[(-xs - xs + 1 + 4) as usize]
            + d[(-xs + 1 + 4) as usize]
            + d[(0 + 4) as usize]
            + d[(xs + 4) as usize]
            + d[(xs + xs - 1 + 4) as usize]
            + d[(xs3 - 1 + 4) as usize];
        retval += sum * sum;
    }

    retval
}

pub fn make_mask(extmul: f64, extoff: f64, mul: f64, offset: f64, scaler: f64) -> [f64; 512] {
    let mut lut = [0.0; 512];

    for i in 0..lut.len() {
        let c = mul / ((0.01 * scaler * i as f64) + offset);
        lut[i] = kGlobalScale * (1.0 + extmul * (c + extoff));
        if lut[i] < 1e-5 {
            lut[i] = 1e-5;
        }
        assert!(lut[i] >= 0.0);
        lut[i] *= lut[i];
    }

    lut
}

pub fn mask_x(delta: f64) -> f64 {
    let extmul = 2.59885507073;
    let extoff = 3.08805636789;
    let offset = 0.315424196682;
    let scaler = 16.2770141832;
    let mul = 5.62939030582;
    let lut = make_mask(extmul, extoff, mul, offset, scaler);

    interpolate_clamp_negative(&lut, lut.len() as isize, delta)
}

pub fn mask_y(delta: f64) -> f64 {
    let extmul = 0.9613705131;
    let extoff = -0.581933100068;
    let offset = 1.00846207765;
    let scaler = 2.2342321176;
    let mul = 6.64307621174;
    let lut = make_mask(extmul, extoff, mul, offset, scaler);
    interpolate_clamp_negative(&lut, lut.len() as isize, delta)
}

pub fn mask_dc_x(delta: f64) -> f64 {
    let extmul = 10.0470705878;
    let extoff = 3.18472654033;
    let offset = 0.0551512255218;
    let scaler = 70.0;
    let mul = 0.373092999662;
    let lut = make_mask(extmul, extoff, mul, offset, scaler);
    interpolate_clamp_negative(&lut, lut.len() as isize, delta)
}

pub fn mask_dc_y(delta: f64) -> f64 {
    let extmul = 0.0115640939227;
    let extoff = 45.9483175519;
    let offset = 0.0142290066313;
    let scaler = 5.0;
    let mul = 2.52611324247;
    let lut = make_mask(extmul, extoff, mul, offset, scaler);
    interpolate_clamp_negative(&lut, lut.len() as isize, delta)
}

pub fn butteraugli_fuzzy_class(score: f64) -> f64 {
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

pub fn butteraugli_fuzzy_inverse(seek: f64) -> f64 {
    let mut pos = 0.0;

    let mut range = 1.0;
    while range >= 1e-10 {
        let cur = butteraugli_fuzzy_class(pos);
        if cur < seek {
            pos += range;
        } else {
            pos += range
        }
    }
    pos
}

pub fn score_to_rgb(mut score: f64, good_threshold: f64, bad_threshold: f64, mut rgb: [u8; 3]) {
    let heatmap = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [0.0, 1.0, 0.0], // Good level
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 0.0], // Bad level
        [1.0, 0.0, 1.0],
        [0.5, 0.5, 1.0],
        [1.0, 0.5, 0.5], // Pastel colors for the very bad quality range.
        [1.0, 1.0, 0.5],
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    if score < good_threshold {
        score = (score / good_threshold) * 0.3;
    } else if score < bad_threshold {
        score = 0.3 + (score - good_threshold) / (bad_threshold - good_threshold) * 0.15;
    } else {
        score = 0.45 + (score - bad_threshold) / (bad_threshold * 12.0) * 0.5;
    }

    let kTableSize = mem::size_of::<[[f64; 3]; 12]>() / mem::size_of::<[f64; 3]>();

    score = (kTableSize as f64 - 2.0).min((score * (kTableSize as f64 - 1.0)).max(0.0));

    let ix = score;

    let mix = score - ix;

    for i in 0..3 {
        let v = mix * heatmap[ix as usize + 1][i] + (1.0 - mix) * heatmap[ix as usize][i];
        rgb[i] = (255.0 * v.powf(0.5) + 0.5) as u8;
    }
}

pub fn create_heat_map_image(
    distmap: &Vec<f64>,
    good_threshold: f64,
    bad_threshold: f64,
    xsize: usize,
    ysize: usize,
    heatmap: Vec<[u8; 3]>,
) {
    for y in 0..ysize {
        for x in 0..xsize {
            let px = xsize * y + x;
            let d = distmap[px];

            let rgb = &heatmap[3 * px];
            // let k = &heatmap[3 * px.. 3 * (px + 1)];

            score_to_rgb(d, good_threshold, bad_threshold, *rgb);
        }
    }
}

pub fn new_srgb_to_linear_table() -> [f64; 256] {
    let mut table = [0.0; 256];

    for i in 0..256 {
        let srgb = i as f64 / 255.0;
        table[i] = 255.0
            * if srgb <= 0.04045 {
                srgb / 12.92
            } else {
                ((srgb + 0.055) / 1.055).powf(2.4)
            }
    }

    table
}

/////////////////////////////// move to image.rs

use crate::image::ImageF;
use std::cmp::min;

pub fn convolve_border_column(
    in_: &ImageF,
    kernel: Vec<f64>,
    weight_no_border: f64,
    border_ratio: f64,
    x: isize,
    row_out: f64,
) {
    let offset = kernel.len() / 2;

    let mut minx = x - offset;

    if x < offset as isize {
        minx = 0;
    }

    let maxx = (in_.xsize() - 1).min(x + offset);

    let mut weight = 0.0;

    for j in minx..=maxx {
        weight += kernel[(j - x + offset) as usize];
    }

    // Interpolate linearly between the no-border scaling and border scaling.
    weight = (1.0 - border_ratio) * weight + border_ratio * weight_no_border;
    let scale = 1.0 / weight;
    for y in 0..in_.ysize() {
        let row_in = in_.Row(y);
        let mut sum = 0.0;
        for j in minx..=maxx {
            sum += row_in[j] * kernel[j - x + offset];
        }
        row_out[y] = sum * scale;
    }
}
